package org.campagnelab.dl.varanalysis.intermediaries;

import org.campagnelab.dl.varanalysis.format.PosRecord;
import org.campagnelab.dl.varanalysis.format.SampleRecord;
import org.campagnelab.dl.varanalysis.storage.AvroVariationParquetReader;
import org.campagnelab.dl.varanalysis.storage.AvroVariationParquetWriter;
import it.unimi.dsi.util.XorShift128PlusRandom;
import org.apache.avro.generic.GenericData;

import java.util.ArrayList;
import java.util.List;
import java.util.Random;

/**
 * Created by rct66 on 5/18/16.
 *
 * the mutator object iterates over a parquet file and creates an additional mutated copy of every record.
 *
 *
 */
public class Mutator extends Intermediary{
    //delta will be halved in homozygous cases (to account for twice the reads at a base)
    //min fraction of bases mutated at a record (ceilinged fraction)
    double deltaSmall = 0.1;
    //max fraction of bases mutated at a record (floored fraction)
    double deltaBig = 1;
    //minimum proportion of counts to presume allele
    double zygHeuristic = 0.1;

    public Mutator(){
    }

    public Mutator(int deltaSmall, int deltaBig, int zygHeuristic){
        this.deltaSmall = deltaSmall;
        this.deltaBig = deltaBig;
        this.zygHeuristic = zygHeuristic;
    }

    public void execute(String in, String out, int blockSize, int pageSize) {
        AvroVariationParquetReader reader = new AvroVariationParquetReader(in);
        AvroVariationParquetWriter writer = new AvroVariationParquetWriter(out,blockSize,pageSize);
        PosRecord pos;
        while(true) {
            //read and write unmutated record
            pos = reader.read();
            if (pos == null)
                break;
            writer.writeRecord(pos);
            //mutate record and write it again
            pos.setMutated(true);
            List<SampleRecord> sampleList = pos.getSamples();
            SampleRecord mutating = sampleList.get(sampleList.size() - 1);
            mutating.setCounts(mutate(mutating.getCounts()));
            pos.setSamples(sampleList);
            //record classes do not deep copy
            writer.writeRecord(pos);
        }
        reader.close();
        writer.close();

    }

    //backward strand appended to forward strand in input and output
    //10 fields: ATCGOATCGO
    //List instead of array because of avro code generation...
    private List<Integer> mutate(List<Integer> counts){
        int[] forward = new int[5];
        int[] backward = new int[5];
        int[] sums = new int[5];

        //fill declared arrays
        int i = 0;
        for (Integer count : counts){
            if (i < 5){
                forward[i] = count;
                sums[i] = count;
            } else {
                backward[i-5] = count;
                sums[i-5] += count;
            }
            i++;
        }
        int maxCount = 0;
        int maxCountIdx = -1;
        int scndCountIdx = -1;
        int numCounts = 0;
        //find highest count idx, second highest count idx, and record number of counts
        for (i = 0; i < 4; i++){
            numCounts += sums[i];
            if (sums[i] > maxCount){
                scndCountIdx = maxCountIdx;
                maxCountIdx = i;
            }
        }
        //no reads whatsoever
        if (maxCountIdx == -1){
            return counts;
        }
        boolean monozygotic;
        //all reads same base, monozygotic
        if (scndCountIdx == -1){
            monozygotic = true;
        } else {
            //see if base with second most reads exceeds heuristic
            monozygotic = (zygHeuristic * numCounts > sums[scndCountIdx]);
        }
        //make rand generator and generate proportion mutating bases
        Random rand = new XorShift128PlusRandom();

        //generate mutation rate
        double delta = deltaSmall + ((deltaBig - deltaSmall) * rand.nextDouble());

        int newBase;
        int oldBase;

        if (monozygotic){
            oldBase = maxCountIdx;
            //generate from non-max bases uniformly
            newBase = rand.nextInt(3);
            if (newBase == oldBase){
                newBase = 3;
            }
            //only one allele mutates, so halve delta when monozygotic
            delta = delta/2;
        } else {
            boolean mutatingAllele = rand.nextBoolean();
            oldBase = mutatingAllele ? maxCountIdx : scndCountIdx;
            newBase = rand.nextInt(3);
            if (newBase == oldBase){
                newBase = 3;
            }
        }
        int oldCount = forward[oldBase];
        for (i = 0; i < oldCount; i++){
            if (rand.nextDouble() < delta) {
                forward[oldBase]--;
                forward[newBase]++;
            }
        }
        oldCount = backward[oldBase];
        for (i = 0; i < oldCount; i++){
            if (rand.nextDouble() < delta) {
                backward[oldBase]--;
                backward[newBase]++;
            }
        }

        //convert two consecutive arrays to a List<Integer>
        List<Integer> newCounts = new ArrayList<Integer>(10);
        for (i = 0; i < 10; i++){
            if (i < 5){
                newCounts.add(forward[i]);
            } else {
                newCounts.add(backward[i-5]);
            }
        }
        return newCounts;
    }

}
