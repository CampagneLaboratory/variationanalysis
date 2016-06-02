package org.campagnelab.dl.varanalysis.intermediaries;

import it.unimi.dsi.logging.ProgressLogger;
import it.unimi.dsi.util.XorShift128PlusRandom;
import org.campagnelab.dl.varanalysis.learning.DetectMutations;
import org.campagnelab.dl.varanalysis.protobuf.BaseInformationRecords;
import org.campagnelab.dl.varanalysis.storage.ParquetPrinter;
import org.campagnelab.dl.varanalysis.storage.RecordReader;
import org.campagnelab.dl.varanalysis.storage.RecordWriter;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.Random;

/**
 * The mutator object iterates over a parquet file and creates an additional mutated copy of every record.
 * <p>
 * Created by rct66 on 5/18/16.
 *
 * @author rct66
 */
public class Mutator extends Intermediary {
    static private Logger LOG = LoggerFactory.getLogger(Mutator.class);

    //delta will be halved in homozygous cases (to account for twice the reads at a base)
    //min fraction of bases mutated at a record (ceilinged fraction)
    double deltaSmall = 0.1;
    //max fraction of bases mutated at a record (floored fraction)
    double deltaBig = 1;
    //minimum proportion of counts to presume allele
    double zygHeuristic = 0.1;
    final String[] STRING = new String[]{"A", "T", "C", "G"};


    public static void main(String[] args) throws IOException{
        //new ParquetPrinter(args[0]).print();
        new Mutator().executeOver(args[0],args[1]);
        new Randomizer().executeOver(args[1],args[2]);
        new ParquetPrinter(args[2]).print();
    }

    public Mutator() {
        setSeed(2323);
    }

    public void setSeed(int seed) {

       rand= new XorShift128PlusRandom(seed);
    }

    public Mutator(int deltaSmall, int deltaBig, int zygHeuristic) {
        this.deltaSmall = deltaSmall;
        this.deltaBig = deltaBig;
        this.zygHeuristic = zygHeuristic;
        setSeed(2323);
    }

    public void execute(String in, String out, int blockSize, int pageSize) throws IOException {
        RecordReader reader = new RecordReader(in);
        RecordWriter writer = new RecordWriter(out, blockSize, pageSize,true);

        //set up logger
        ProgressLogger pgReadWrite = new ProgressLogger(LOG);
        pgReadWrite.itemsName = "mutation";
        pgReadWrite.expectedUpdates = reader.getTotalRecords();
        pgReadWrite.displayFreeMemory = true;
        pgReadWrite.start();

        for (BaseInformationRecords.BaseInformation base : reader) {
            BaseInformationRecords.BaseInformation.Builder baseBuilder = base.toBuilder();

            //set is tumor flags (false sample 0 (germline), true sample 1 (somatic)
            baseBuilder.setSamples(0,baseBuilder.getSamples(0).toBuilder().setIsTumor(false));
            baseBuilder.setSamples(1,baseBuilder.getSamples(1).toBuilder().setIsTumor(true));
            writer.writeRecord(baseBuilder.build());

            //mutate record and write it again
            writer.writeRecord(mutate(baseBuilder));
            pgReadWrite.update();
        }
        pgReadWrite.stop();
        reader.close();
        writer.close();

    }

    Random rand;


    //backward strand appended to forward strand in input and output
    //10 fields: ATCGOATCGO
    //List instead of array because of avro code generation...
    protected BaseInformationRecords.BaseInformation mutate(BaseInformationRecords.BaseInformation.Builder baseBuild) {
        baseBuild.setMutated(true);
        int[] forward = new int[5];
        int[] backward = new int[5];
        int[] sums = new int[5];
        BaseInformationRecords.SampleInfo somatic = baseBuild.getSamples(1);
        //fill declared arrays
        int i = 0;
        for (BaseInformationRecords.CountInfo count : somatic.getCountsList()) {
            forward[i] = count.getGenotypeCountForwardStrand();
            backward[i] = count.getGenotypeCountReverseStrand();
            sums[i] = forward[i] + backward[i];
            i++;
        }
        int maxCount = 0;
        int maxCountIdx = -1;
        int scndCountIdx = -1;
        int numCounts = 0;
        //find highest count idx, second highest count idx, and record number of counts
        for (i = 0; i < 4; i++) {
            numCounts += sums[i];
            if (sums[i] > maxCount) {
                scndCountIdx = maxCountIdx;
                maxCountIdx = i;
                maxCount = sums[i];
            }
        }
        //no reads whatsoever
        if (maxCountIdx == -1) {
            return baseBuild.build();
        }
        boolean monozygotic;
        //all reads same base, monozygotic
        if (scndCountIdx == -1) {
            monozygotic = true;
        } else {
            //see if base with second most reads exceeds heuristic
            monozygotic = (zygHeuristic * numCounts > sums[scndCountIdx]);
        }
        //make rand generator and generate proportion mutating bases
        //generate mutation rate
        double delta = deltaSmall + ((deltaBig - deltaSmall) * rand.nextDouble());
        double deltaOrig = delta;

        int newBase;
        int oldBase;

        if (monozygotic) {
            oldBase = maxCountIdx;
            //generate from non-max bases uniformly
            newBase = rand.nextInt(3);
            if (newBase == oldBase) {
                newBase = 3;
            }
            //only one allele mutates, so halve delta when monozygotic
            delta = delta / 2;
        } else {
            boolean mutatingAllele = rand.nextBoolean();
            oldBase = mutatingAllele ? maxCountIdx : scndCountIdx;
            newBase = rand.nextInt(3);
            if (newBase == oldBase) {
                newBase = 3;
            }
        }
        int oldCount = forward[oldBase];
        for (i = 0; i < oldCount; i++) {
            if (rand.nextDouble() < delta) {
                forward[oldBase]--;
                forward[newBase]++;
            }
        }
        oldCount = backward[oldBase];
        for (i = 0; i < oldCount; i++) {
            if (rand.nextDouble() < delta) {
                backward[oldBase]--;
                backward[newBase]++;
            }
        }
        //write to respective builders and return rebuild
        BaseInformationRecords.SampleInfo.Builder somaticBuild = somatic.toBuilder();
        i = 0;
        for (BaseInformationRecords.CountInfo count : somaticBuild.getCountsList()) {
            BaseInformationRecords.CountInfo.Builder countBuild = count.toBuilder();
            countBuild.setGenotypeCountForwardStrand(forward[i]);
            countBuild.setGenotypeCountReverseStrand(backward[i]);
            somaticBuild.setCounts(i, countBuild);
            i++;
        }
        baseBuild.setSamples(1, somaticBuild);
        baseBuild.setFrequencyOfMutation((float) deltaOrig);
        String newBaseString = STRING[newBase];
        baseBuild.setMutatedBase(newBaseString);
        baseBuild.setIndexOfMutatedBase(newBase);
        return baseBuild.build();
    }

}