package org.campagnelab.dl.somatic.intermediaries;

import it.unimi.dsi.fastutil.ints.IntArrayList;
import it.unimi.dsi.fastutil.ints.IntArraySet;
import it.unimi.dsi.util.XorShift1024StarRandom;
import org.campagnelab.dl.varanalysis.protobuf.BaseInformationRecords;

import java.util.Collections;

/**
 * Created by fac2003 on 7/19/16.
 */
public class TwoSampleCanonicalSimulationStrategy implements SimulationStrategy {
    private XorShift1024StarRandom randomGenerator;
    private double canonThreshold;

    FirstSimulationStrategy firstSimulationStrategy;

    IntArrayList genotypeCounts0 = new IntArrayList();
    IntArrayList genotypeCounts1 = new IntArrayList();
    IntArrayList sortingPermutationGenotypeCounts0 = new IntArrayList();
    IntArrayList sortingPermutationGenotypeCounts1 = new IntArrayList();

    public TwoSampleCanonicalSimulationStrategy() {


    }

    public TwoSampleCanonicalSimulationStrategy(double deltaSmall, double deltaBig, double zygHeuristic, long seed, double canonThreshold) {
        setup(deltaSmall, deltaBig, zygHeuristic, seed, canonThreshold);
    }

    public void setup(double deltaSmall, double deltaBig, double zygHeuristic, long seed, double canonThreshold) {
        this.canonThreshold = canonThreshold;
        firstSimulationStrategy = new FirstSimulationStrategy();
        firstSimulationStrategy.setup(deltaSmall, deltaBig, zygHeuristic, seed, 0);
        setSeed(seed);
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
        int numAlleles0 = sumGenotype90P(sumCount0, genotypeCounts0, sortingPermutationGenotypeCounts0);
        int numAlleles1 = sumGenotype90P(sumCount1, genotypeCounts1, sortingPermutationGenotypeCounts1);
        if (numAlleles0 != numAlleles1 || numAlleles0 > 2) {
            // we won't make this somatic because the site is not like the germline site we expect.
            // by not making a somatic variant we will prevent the model learning that such sites are valid predictions.
            makeSomatic = false;
        }
        //find unique allele indices to make sure there aren't too many
        IntArraySet alleleIndices = new IntArraySet();
        alleleIndices.add(sortingPermutationGenotypeCounts0.getInt(0));
        alleleIndices.add(sortingPermutationGenotypeCounts1.getInt(0));
        //heterozygous case
        if (numAlleles0 > 1) {
            alleleIndices.add(sortingPermutationGenotypeCounts0.getInt(1));
            alleleIndices.add(sortingPermutationGenotypeCounts1.getInt(1));

        }
        if (alleleIndices.size() > numAlleles0) {
            // this site is not canonical.
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