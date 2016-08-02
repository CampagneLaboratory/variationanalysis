package org.campagnelab.dl.varanalysis.intermediaries;

import it.unimi.dsi.fastutil.ints.IntArrayList;
import it.unimi.dsi.util.XoRoShiRo128PlusRandom;
import org.campagnelab.dl.varanalysis.protobuf.BaseInformationRecords;

import java.util.Collections;
import java.util.Date;

/**
 * Created by fac2003 on 7/19/16.
 */
public class SimulationStrategyImplTrio implements SimulationStrategy {
    private XoRoShiRo128PlusRandom randomGenerator;
    private long seed;

    public SimulationStrategyImplTrio(long seed) {
        this.seed = seed;
    }

    public SimulationStrategyImplTrio() {
        this.seed = new Date().getTime();
        randomGenerator = new XoRoShiRo128PlusRandom(seed);
        firstSimulationStrategy = new FirstSimulationStrategy(seed);
    }

    FirstSimulationStrategy firstSimulationStrategy;

    IntArrayList genotypeCounts0 = new IntArrayList();
    IntArrayList genotypeCounts1 = new IntArrayList();
    IntArrayList genotypeCounts2 = new IntArrayList();
    IntArrayList sortingPermutationGenotypeCounts0 = new IntArrayList();
    IntArrayList sortingPermutationGenotypeCounts1 = new IntArrayList();
    IntArrayList sortingPermutationGenotypeCounts2 = new IntArrayList();

    @Override
    public BaseInformationRecords.BaseInformation mutate(boolean makeSomatic,
                                                         BaseInformationRecords.BaseInformation record,
                                                         BaseInformationRecords.SampleInfo unused1,
                                                         BaseInformationRecords.SampleInfo unused2,
                                                         SimulationCharacteristics sim) {
        //overwrite parameters for trio, this is workaround so that strategy interface can remain unchanged
        BaseInformationRecords.SampleInfo germlineSample = record.getSamples(2);
        BaseInformationRecords.SampleInfo fatherSample = record.getSamples(0);
        BaseInformationRecords.SampleInfo motherSample = record.getSamples(1);

        //note that from here, 0=germ, 1=father, 2=mother
        int sumCountGerm = getSumCounts(germlineSample, genotypeCounts0);
        int sumCountFather = getSumCounts(fatherSample, genotypeCounts1);
        int sumCountMother = getSumCounts(motherSample, genotypeCounts2);
        prepareSorted(genotypeCounts0, sortingPermutationGenotypeCounts0);
        prepareSorted(genotypeCounts1, sortingPermutationGenotypeCounts1);
        prepareSorted(genotypeCounts2, sortingPermutationGenotypeCounts2);

        int numAlleles0 = sumGenotype90P(sumCountGerm, genotypeCounts0, sortingPermutationGenotypeCounts0);
        int numAlleles1 = sumGenotype90P(sumCountFather, genotypeCounts1, sortingPermutationGenotypeCounts1);
        int numAlleles2 = sumGenotype90P(sumCountMother, genotypeCounts1, sortingPermutationGenotypeCounts2);

        // check that alleles can possibly come from parents

        //define genotypes, eg AA or AB.
        int child1 = sortingPermutationGenotypeCounts0.getInt(0);
        int child2 = (numAlleles0>1)?sortingPermutationGenotypeCounts0.getInt(1):child1;
        int father1 = sortingPermutationGenotypeCounts1.getInt(0);
        int father2 = (numAlleles1>1)?sortingPermutationGenotypeCounts1.getInt(1):father1;
        int mother1 = sortingPermutationGenotypeCounts2.getInt(0);;
        int mother2 = (numAlleles2>1)?sortingPermutationGenotypeCounts2.getInt(1):mother1;

        //first, check that germline/child doesn't have too many alleles (or is not designated for mutation). if not, then
        //determine if the child's genotype is possible, allowing either first allele from father or from mother
        if (numAlleles0 > 2 || (!makeSomatic)){
            makeSomatic = false;
        } else {
            makeSomatic = isMendelian(child1,child2,father1,father2,mother1,mother2);
        }
        return firstSimulationStrategy.mutate(makeSomatic, record, null, null, null);

    }
    public static boolean isMendelian(int child1, int child2, int father1, int father2, int mother1, int mother2){
        boolean firstFromFather = ((child1 == father1) || (child1 == father2)) && ((child2 == mother1) || (child2 == mother2));
        boolean firstFromMother = ((child1 == mother1) || (child1 == mother2)) && ((child2 == father1) || (child2 == father2));
        return (firstFromFather ||  firstFromMother);
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
