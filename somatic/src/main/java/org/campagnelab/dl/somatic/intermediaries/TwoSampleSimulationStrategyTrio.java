package org.campagnelab.dl.somatic.intermediaries;

import it.unimi.dsi.fastutil.ints.IntArrayList;
import it.unimi.dsi.util.XorShift1024StarRandom;
import org.campagnelab.dl.varanalysis.protobuf.BaseInformationRecords;

import java.util.Collections;

/**
 * Created by rct66 on 12/12/16.
 * Although called onesample after the similar strategy used for pair models, this stragey checks that just the mother
 * and father samples (not the patient sample) are canonical before mutating. Mendelian rules do not have to be followed.
 */
public class TwoSampleSimulationStrategyTrio implements SimulationStrategy {
    private XorShift1024StarRandom randomGenerator;
    private double canonThreshold;

    public TwoSampleSimulationStrategyTrio(long seed) {
        setSeed(seed);
    }

    public TwoSampleSimulationStrategyTrio() {

    }

    @Override
    public void setup(double deltaSmall, double deltaBig, double zygHeuristic, long seed, double canonThreshold) {
        setSeed(seed);
        this.canonThreshold = canonThreshold;
        firstSimulationStrategy = new FirstSimulationStrategy();
        firstSimulationStrategy.setup(deltaSmall, deltaBig, zygHeuristic, seed,0);
    }

    FirstSimulationStrategy firstSimulationStrategy;

    IntArrayList genotypeCounts0 = new IntArrayList();
    IntArrayList genotypeCounts1 = new IntArrayList();
    IntArrayList sortingPermutationGenotypeCounts0 = new IntArrayList();
    IntArrayList sortingPermutationGenotypeCounts1 = new IntArrayList();

    @Override
    public int numberOfSamplesSupported() {
        return 3;
    }

    @Override
    public BaseInformationRecords.BaseInformation mutate(boolean makeSomatic,
                                                         BaseInformationRecords.BaseInformation record,
                                                         BaseInformationRecords.SampleInfo unused1,
                                                         BaseInformationRecords.SampleInfo unused2,
                                                         SimulationCharacteristics sim) {
        //overwrite parameters for trio, this is workaround so that strategy interface can remain unchanged
        BaseInformationRecords.SampleInfo fatherSample = record.getSamples(0);
        BaseInformationRecords.SampleInfo motherSample = record.getSamples(1);

        //note that from here, 0=germ, 1=father, 2=mother
        int sumCountFather = getSumCounts(fatherSample, genotypeCounts0);
        int sumCountMother = getSumCounts(motherSample, genotypeCounts1);
        prepareSorted(genotypeCounts0, sortingPermutationGenotypeCounts0);
        prepareSorted(genotypeCounts1, sortingPermutationGenotypeCounts1);

        int numAlleles0 = sumGenotype90P(sumCountFather, genotypeCounts0, sortingPermutationGenotypeCounts0);
        int numAlleles1 = sumGenotype90P(sumCountMother, genotypeCounts1, sortingPermutationGenotypeCounts1);

        if (numAlleles0 > 2 || numAlleles1 > 2) {

            // we won't make this somatic because the either the father or mother sample sample is not diploid.
            // by not making a somatic variant we will prevent the model learning that such sites are valid predictions.
            makeSomatic = false;
        }

        return firstSimulationStrategy.mutate(makeSomatic, record, null, null, sim);

    }



    @Override
    public void setSeed(long seed) {
        randomGenerator = new XorShift1024StarRandom(seed);
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