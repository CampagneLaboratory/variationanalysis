package org.campagnelab.dl.somatic.intermediaries;

import it.unimi.dsi.fastutil.ints.IntArrayList;
import it.unimi.dsi.util.XorShift1024StarRandom;
import org.campagnelab.dl.varanalysis.protobuf.BaseInformationRecords;

import java.util.Collections;

/**
 * A strategy that only considers the germline sample to determine if the site is canonical.
 * Created by fac2003 on 7/19/16.
 */
public class OneSampleCanonicalSimulationStrategy implements SimulationStrategy {
    private XorShift1024StarRandom randomGenerator;
    private double canonThreshold;

    FirstSimulationStrategy firstSimulationStrategy;

    IntArrayList genotypeCounts0 = new IntArrayList();
    IntArrayList genotypeCounts1 = new IntArrayList();
    IntArrayList sortingPermutationGenotypeCounts0 = new IntArrayList();
    IntArrayList sortingPermutationGenotypeCounts1 = new IntArrayList();

    public OneSampleCanonicalSimulationStrategy() {

    }

    @Override
    public void setup(double deltaSmall, double deltaBig, double zygHeuristic, long seed, double canonThreshold) {
        firstSimulationStrategy = new FirstSimulationStrategy();
        firstSimulationStrategy.setup(deltaSmall, deltaBig, zygHeuristic, seed, 0);
        this.canonThreshold = canonThreshold;
        setSeed(seed);
    }

    @Override
    public int numberOfSamplesSupported() {
        return 2;
    }

    @Override
    public BaseInformationRecords.BaseInformation mutate(boolean makeSomatic,
                                                         BaseInformationRecords.BaseInformation record,
                                                         BaseInformationRecords.SampleInfo germlineSample,
                                                         BaseInformationRecords.SampleInfo otherSample,
                                                         SimulationCharacteristics sim) {

        int sumCount0 = getSumCounts(germlineSample, genotypeCounts0);
        int sumCount1 = getSumCounts(otherSample, genotypeCounts1);
        prepareSorted(genotypeCounts0, sortingPermutationGenotypeCounts0);
        prepareSorted(genotypeCounts1, sortingPermutationGenotypeCounts1);
        int numAllelesGermline = sumGenotype90P(sumCount0, genotypeCounts0, sortingPermutationGenotypeCounts0);

        if (numAllelesGermline > 2) {

            // we won't make this somatic because the the germline sample is not diploid.
            // by not making a somatic variant we will prevent the model learning that such sites are valid predictions.
            makeSomatic = false;
        }

        return firstSimulationStrategy.mutate(makeSomatic, record, germlineSample, otherSample, null);
    }


    @Override
    public void setSeed(long seed) {
        randomGenerator = new XorShift1024StarRandom(seed);
        firstSimulationStrategy.setSeed(seed);
    }

    private void prepareSorted(IntArrayList original, IntArrayList permutation) {

        permutation.clear();
        for (int i = 0; i < original.size(); i++) {
            permutation.add(i);
        }
        // sort using the counts of the original
        Collections.sort(permutation, (o1, o2) -> original.getInt(o2) - original.getInt(o1));

    }

    /**
     * Calculate the number of genotypes with the most counts that together achieve 90% of total counts.
     *
     * @param sumCount
     * @param genotypeCounts                   the counts array, in genotype order
     * @param sortingPermutationGenotypeCounts the permutation array, from cumulative order to genotype order
     * @return
     */
    private int sumGenotype90P(int sumCount, IntArrayList genotypeCounts, IntArrayList sortingPermutationGenotypeCounts) {
        int cumulative = 0;
        int index = 0;
        for (int i = 0; i < sortingPermutationGenotypeCounts.size(); i++) {
            int count = getSortedCountAtIndex(i, genotypeCounts, sortingPermutationGenotypeCounts);
            cumulative += count;
            index++;
            if (cumulative > (sumCount * canonThreshold)) return index;

        }
        return index;
    }

    private int getSortedCountAtIndex(int i, IntArrayList genotypeCounts, IntArrayList sortingPermutationGenotypeCounts) {
        return genotypeCounts.getInt(sortingPermutationGenotypeCounts.getInt(i));
    }

    private int getSumCounts(BaseInformationRecords.SampleInfo sample, IntArrayList genotypeCounts) {
        genotypeCounts.clear();
        int sumCounts = 0;
        for (BaseInformationRecords.CountInfo c : sample.getCountsList()) {
            final int count = c.getGenotypeCountForwardStrand() + c.getGenotypeCountReverseStrand();
            genotypeCounts.add(count);
            sumCounts += count;
        }
        return sumCounts;
    }


}