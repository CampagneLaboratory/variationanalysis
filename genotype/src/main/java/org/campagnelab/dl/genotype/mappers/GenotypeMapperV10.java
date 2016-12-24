package org.campagnelab.dl.genotype.mappers;

/**
 * With combined layer. Based on V9.
 */
public class GenotypeMapperV10 extends GenotypeMapperV9 {
    public GenotypeMapperV10() {
        super();
        sortCounts = true;
        withDistinctAlleleCounts = false;
        withCombinedLayer = true;
    }

}
