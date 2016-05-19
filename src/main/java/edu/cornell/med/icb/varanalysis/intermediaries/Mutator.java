package edu.cornell.med.icb.varanalysis.intermediaries;

import edu.cornell.med.icb.varanalysis.format.PosRecord;
import edu.cornell.med.icb.varanalysis.format.SampleRecord;
import edu.cornell.med.icb.varanalysis.storage.AvroVariationParquetReader;
import edu.cornell.med.icb.varanalysis.storage.AvroVariationParquetWriter;
import it.unimi.dsi.util.XorShift128PlusRandom;

import java.util.List;
import java.util.Random;

/**
 * Created by rct66 on 5/18/16.
 *
 * the mutator object iterates over a parquet file and creates an additional mutated copy of every record.
 * Mutations
 *
 */
public class Mutator{
    //delta will be halved in homozygous cases (to account for twice the reads at a base)
    //min fraction of bases mutated at a record (ceilinged fraction)
    double deltaSmall = 0.1;
    //max fraction of bases mutated at a record (floored fraction)
    double deltaBig = 1;
    //minimum proportion of counts to presume allele
    double zygHeuristic = 0.1;



    void process(String in, String out, int blockSize, int pageSize) {
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
                sums[i] += count;
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



        //monozygotic case first:
        if (monozygotic){
            //generate from non-max bases uniformly
            int newBase = rand.nextInt(3);
            if (newBase == maxCountIdx){
                newBase = 3;
            }
            int add = 0;
            int sub = 0;
            for (i = 0; i < forward[newBase]; i++){
                //if rand.next
            }
        }

        return null;
    }

}
