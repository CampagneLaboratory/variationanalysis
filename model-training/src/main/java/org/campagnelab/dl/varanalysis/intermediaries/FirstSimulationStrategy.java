package org.campagnelab.dl.varanalysis.intermediaries;

import it.unimi.dsi.fastutil.ints.IntArrayList;
import it.unimi.dsi.fastutil.objects.ObjectArrayList;
import it.unimi.dsi.util.XorShift1024StarRandom;
import org.campagnelab.dl.model.utils.ProtoPredictor;
import org.campagnelab.dl.varanalysis.protobuf.BaseInformationRecords;

import java.util.Collections;
import java.util.Date;
import java.util.List;
import java.util.Random;

/**
 * Introduce a strategy to introduce mutations.
 * Created by fac2003 on 7/19/16.
 */
public class FirstSimulationStrategy implements SimulationStrategy {
    long seed;

    public FirstSimulationStrategy(long seed) {
        this.seed = seed;

        rand = new XorShift1024StarRandom(seed);
    }

    public FirstSimulationStrategy(int deltaSmall, int deltaBig, int zygHeuristic, long seed) {
        this(new Date().getTime());

        this.deltaSmall = deltaSmall;
        this.deltaBig = deltaBig;
        this.zygHeuristic = zygHeuristic;
    }


    //delta will be halved in homozygous cases (to account for twice the reads at a base)
    //min fraction of bases mutated at a record (ceilinged fraction)
    double deltaSmall = 0.0;
    //max fraction of bases mutated at a record (floored fraction)
    double deltaBig = 1;
    //minimum proportion of counts to presume allele
    double zygHeuristic = 0.1;
    final String[] STRING = new String[]{"A", "T", "C", "G"};
    Random rand;

    @Override
    public BaseInformationRecords.BaseInformation mutate(boolean makeSomatic,
                                                         BaseInformationRecords.BaseInformation record,
                                                         BaseInformationRecords.SampleInfo germlineSample,
                                                         BaseInformationRecords.SampleInfo otherSample,
                                                         SimulationCharacteristics sim) {

        BaseInformationRecords.BaseInformation.Builder baseBuild = record.toBuilder();
        baseBuild.setMutated(makeSomatic);
        int numSamples = record.getSamplesCount();
        for (int i = 0; i < numSamples-1; i++){
            baseBuild.setSamples(i, record.getSamples(i).toBuilder().setIsTumor(false));
        }
        baseBuild.setSamples(numSamples-1, record.getSamples(numSamples-1).toBuilder().setIsTumor(true));
        if (!makeSomatic) {
            // don't change counts if we are not to make the second sample somatic.
            return baseBuild.build();
        }
        BaseInformationRecords.SampleInfo somatic = baseBuild.getSamples(numSamples-1);
        int numGenos = somatic.getCountsList().size();
        int[] forward = new int[numGenos];
        int[] backward = new int[numGenos];
        int[] sums = new int[numGenos];
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
        for (i = 0; i < numGenos; i++) {
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
        int otherBase = -1;
        if (monozygotic) {
            oldBase = maxCountIdx;
            //generate from non-max bases uniformly

            //only one allele mutates, so halve delta when monozygotic
            delta = delta / 2;
            newBase = rand.nextInt(numGenos - 2);
        } else {
            boolean mutatingAllele = rand.nextBoolean();
            oldBase = mutatingAllele ? maxCountIdx : scndCountIdx;
            otherBase = !mutatingAllele ? maxCountIdx : scndCountIdx;
            newBase = rand.nextInt(numGenos - 3);

        }

        if (newBase == oldBase) {
            //replace self case
            newBase = numGenos - 2;
        } else if (newBase == otherBase){
            //replace
            newBase = numGenos - 3;
        } else if (newBase == 4) {
            //replace genotype N case
            newBase = numGenos - 1;
        }
        int fMutCount = 0;
        int oldCount = forward[oldBase];
        for (i = 0; i < oldCount; i++) {
            if (rand.nextDouble() < delta) {
                forward[oldBase]--;
                forward[newBase]++;
                fMutCount++;

            }
        }
        int bMutCount = 0;
        oldCount = backward[oldBase];
        for (i = 0; i < oldCount; i++) {
            if (rand.nextDouble() < delta) {
                backward[oldBase]--;
                backward[newBase]++;
                bMutCount++;
            }
        }
        //write to respective builders and return rebuild
        BaseInformationRecords.SampleInfo.Builder somaticBuild = somatic.toBuilder();



        //generate mutated numVariationsLists score lists (some boilerplate here...)
        //will add 1 to mutated variant counts using mutateIntegerListsVarAdd
        List<Integer> fromVC = new ObjectArrayList<Integer>();
        List<Integer> toVC = new ObjectArrayList<Integer>();
        if (somaticBuild.getCounts(oldBase).getNumVariationsInReadsCount() > 0) {

            //forward strand
            fromVC.addAll(ProtoPredictor.expandFreq(somaticBuild.getCounts(oldBase).getNumVariationsInReadsList()));
            toVC.addAll(ProtoPredictor.expandFreq(somaticBuild.getCounts(newBase).getNumVariationsInReadsList()));
            mutateIntegerListsVarAdd(fMutCount, fromVC, toVC,somaticBuild.getCounts(oldBase).getMatchesReference(),somaticBuild.getCounts(newBase).getMatchesReference());
        }

        //generate mutated insert sizes lists (some boilerplate here...)
        List<Integer> fromIS = new ObjectArrayList<Integer>();
        List<Integer> toIS = new ObjectArrayList<Integer>();
        if (somaticBuild.getCounts(oldBase).getNumVariationsInReadsCount() > 0) {

            //forward strand
            fromIS.addAll(ProtoPredictor.expandFreq(somaticBuild.getCounts(oldBase).getInsertSizesList()));
            toIS.addAll(ProtoPredictor.expandFreq(somaticBuild.getCounts(newBase).getInsertSizesList()));
            mutateIntegerLists(fMutCount, fromIS, toIS);
        }


        //generate mutated quality score lists (some boilerplate here...)
        //get old list of from scores
        List<Integer> fromForward = new ObjectArrayList<Integer>();
        List<Integer> fromBackward = new ObjectArrayList<Integer>();
        List<Integer> toForward = new ObjectArrayList<Integer>();
        List<Integer> toBackward = new ObjectArrayList<Integer>();
        if (somaticBuild.getCounts(oldBase).getQualityScoresForwardStrandCount() > 0 && somaticBuild.getCounts(oldBase).getQualityScoresReverseStrandCount() > 0) {

            //forward strand
            fromForward.addAll(ProtoPredictor.expandFreq(somaticBuild.getCounts(oldBase).getQualityScoresForwardStrandList()));
            toForward.addAll(ProtoPredictor.expandFreq(somaticBuild.getCounts(newBase).getQualityScoresForwardStrandList()));
            mutateIntegerLists(fMutCount, fromForward, toForward);

            //reverse strand
            fromBackward.addAll(ProtoPredictor.expandFreq(somaticBuild.getCounts(oldBase).getQualityScoresReverseStrandList()));
            toBackward.addAll(ProtoPredictor.expandFreq(somaticBuild.getCounts(newBase).getQualityScoresReverseStrandList()));
            mutateIntegerLists(bMutCount, fromBackward, toBackward);

        }

        //generate mutated readIndex lists
        List<Integer> fromForwardR = new ObjectArrayList<Integer>();
        List<Integer> fromBackwardR = new ObjectArrayList<Integer>();
        List<Integer> toForwardR = new ObjectArrayList<Integer>();
        List<Integer> toBackwardR = new ObjectArrayList<Integer>();
        if (somaticBuild.getCounts(oldBase).getReadIndicesForwardStrandCount() > 0 && somaticBuild.getCounts(oldBase).getReadIndicesReverseStrandCount() > 0) {

            //forward strand
            fromForwardR.addAll(ProtoPredictor.expandFreq(somaticBuild.getCounts(oldBase).getReadIndicesForwardStrandList()));
            toForwardR.addAll(ProtoPredictor.expandFreq(somaticBuild.getCounts(newBase).getReadIndicesForwardStrandList()));
            mutateIntegerLists(fMutCount, fromForwardR, toForwardR);

            //reverse strand
            fromBackwardR.addAll(ProtoPredictor.expandFreq(somaticBuild.getCounts(oldBase).getReadIndicesReverseStrandList()));
            toBackwardR.addAll(ProtoPredictor.expandFreq(somaticBuild.getCounts(newBase).getReadIndicesReverseStrandList()));
            mutateIntegerLists(bMutCount, fromBackwardR, toBackwardR);


        }


        i = 0;
        for (BaseInformationRecords.CountInfo count : somaticBuild.getCountsList()) {
            BaseInformationRecords.CountInfo.Builder countBuild = count.toBuilder();
            countBuild.setGenotypeCountForwardStrand(forward[i]);
            countBuild.setGenotypeCountReverseStrand(backward[i]);
            if (i == oldBase) {
                //replace quality scores
                countBuild.clearQualityScoresForwardStrand();
                countBuild.clearQualityScoresReverseStrand();
                countBuild.addAllQualityScoresForwardStrand(ProtoPredictor.compressFreq(fromForward));
                countBuild.addAllQualityScoresReverseStrand(ProtoPredictor.compressFreq(fromBackward));

                //replace readIndices
                countBuild.clearReadIndicesForwardStrand();
                countBuild.clearReadIndicesReverseStrand();
                countBuild.addAllReadIndicesForwardStrand(ProtoPredictor.compressFreq(fromForwardR));
                countBuild.addAllReadIndicesReverseStrand(ProtoPredictor.compressFreq(fromBackwardR));

                //replace numVars
                countBuild.clearNumVariationsInReads();
                countBuild.addAllNumVariationsInReads(ProtoPredictor.compressFreq(fromVC));

                //replace insert sizes
                countBuild.clearInsertSizes();
                countBuild.addAllInsertSizes(ProtoPredictor.compressFreq(fromIS));


            } else if (i == newBase) {
                //replace quality scores
                countBuild.clearQualityScoresForwardStrand();
                countBuild.clearQualityScoresReverseStrand();
                countBuild.addAllQualityScoresForwardStrand(ProtoPredictor.compressFreq(toForward));
                countBuild.addAllQualityScoresReverseStrand(ProtoPredictor.compressFreq(toBackward));

                //replace readIndices
                countBuild.clearReadIndicesForwardStrand();
                countBuild.clearReadIndicesReverseStrand();
                countBuild.addAllReadIndicesForwardStrand(ProtoPredictor.compressFreq(toForwardR));
                countBuild.addAllReadIndicesReverseStrand(ProtoPredictor.compressFreq(toBackwardR));
                baseBuild.setMutatedBase(count.getToSequence());

                //replace numVars
                countBuild.clearNumVariationsInReads();
                countBuild.addAllNumVariationsInReads(ProtoPredictor.compressFreq(toVC));

                //replace insert sizes
                countBuild.clearInsertSizes();
                countBuild.addAllInsertSizes(ProtoPredictor.compressFreq(toIS));


            }
            somaticBuild.setCounts(i, countBuild);
            i++;
        }
        baseBuild.setSamples(numSamples-1, somaticBuild);
        baseBuild.setFrequencyOfMutation((float) deltaOrig);
        // String newBaseString = newBase<STRING.length? STRING[newBase]:"N";
        //baseBuild.setMutatedBase(newBaseString);
        baseBuild.setIndexOfMutatedBase(newBase);
        return baseBuild.build();
    }

    @Override
    public void setSeed(int seed) {
        rand = new XorShift1024StarRandom(seed);
    }

    private void mutateIntegerLists(int fMutCount, List<Integer> source, List<Integer> dest) {
        Collections.shuffle(source, rand);
        dest.addAll(source.subList(0, fMutCount));
        List<Integer> tmp = new IntArrayList(source.subList(fMutCount, source.size()));
        source.clear();
        source.addAll(tmp);

    }


    private void mutateIntegerListsVarAdd(int fMutCount, List<Integer> source, List<Integer> dest, boolean sourceIsRef, boolean destIsRef) {
        //compute increment (-1,0,or 1)
        int increment = 0;
        if (sourceIsRef) increment++;
        if (destIsRef) increment--;

        int newDestIndex = dest.size();
        Collections.shuffle(source, rand);
        dest.addAll(source.subList(0, fMutCount));
        for (int i = newDestIndex; i < dest.size(); i++){
            dest.set(i,dest.get(i)+increment);
        }
        List<Integer> tmp = new IntArrayList(source.subList(fMutCount, source.size()));
        source.clear();
        source.addAll(tmp);

    }

}



