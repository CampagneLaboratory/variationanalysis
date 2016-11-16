package org.campagnelab.dl.varanalysis.intermediaries;

import it.unimi.dsi.fastutil.ints.IntArrayList;
import it.unimi.dsi.util.XoRoShiRo128PlusRandom;
import org.campagnelab.dl.varanalysis.protobuf.BaseInformationRecords;

import java.util.Collections;
import java.util.Date;

/**
 * Created by fac2003 on 7/19/16.
 */
public class SimulationStrategyImpl implements SimulationStrategy {
    private XoRoShiRo128PlusRandom randomGenerator;
    private long seed;

    public SimulationStrategyImpl(long seed) {
        this.seed = seed;
    }

    public SimulationStrategyImpl() {
        this.seed = new Date().getTime();
        randomGenerator = new XoRoShiRo128PlusRandom(seed);
        firstSimulationStrategy = new FirstSimulationStrategy(seed);
    }

    FirstSimulationStrategy firstSimulationStrategy;

    IntArrayList genotypeCounts0 = new IntArrayList();
    IntArrayList genotypeCounts1 = new IntArrayList();
    IntArrayList sortingPermutationGenotypeCounts0 = new IntArrayList();
    IntArrayList sortingPermutationGenotypeCounts1 = new IntArrayList();

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
        // check that alleles are the same genotype index, or that we have a tie in counts:
        for (int i = 0; i < numAlleles0; i++) {
            if (sortingPermutationGenotypeCounts0.getInt(i) != sortingPermutationGenotypeCounts1.getInt(i) &&
                    getSortedCountAtIndex(i, genotypeCounts0, sortingPermutationGenotypeCounts0) != getSortedCountAtIndex(i, genotypeCounts1, sortingPermutationGenotypeCounts1)) {
                // same number of alleles, but not the same genotypes.
                makeSomatic = false;
            }
        }
        return firstSimulationStrategy.mutate(makeSomatic, record, germlineSample, otherSample, null);
    }

    @Override
    public void setSeed(int seed) {
        firstSimulationStrategy = new FirstSimulationStrategy(seed);
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
     * @param sortingPermutationGenotypeCounts the permutation array, from cumlative order to genotype order
     * @return
     */
    private int sumGenotype90P(int sumCount, IntArrayList genotypeCounts, IntArrayList sortingPermutationGenotypeCounts) {
        int cumulative = 0;
        int index = 0;
        for (int i = 0; i < sortingPermutationGenotypeCounts.size(); i++) {
            int count = getSortedCountAtIndex(i, genotypeCounts, sortingPermutationGenotypeCounts);
            cumulative += count;
            index++;
            if (cumulative > (sumCount * 9 / 10)) return index;

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
