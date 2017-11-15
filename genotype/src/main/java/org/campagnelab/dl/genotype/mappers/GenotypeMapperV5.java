package org.campagnelab.dl.genotype.mappers;

/**
 * This mapper sorts counts and predicts DistinctAlleleCounts.
 */
public class GenotypeMapperV5 extends GenotypeMapperV4 {
    public GenotypeMapperV5() {
        withDistinctAlleleCounts=true;
    }
}
