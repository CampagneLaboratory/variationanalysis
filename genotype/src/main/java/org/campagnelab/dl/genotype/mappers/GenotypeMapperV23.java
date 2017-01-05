package org.campagnelab.dl.genotype.mappers;

/**
 * V122 combined with hasIsVariantLabelMapper.
 */
public class GenotypeMapperV23 extends GenotypeMapperV22 {


    public GenotypeMapperV23() {
        super();
        sortCounts = true;
        withDistinctAlleleCounts = false;
        withCombinedLayer = true;
        hasIsVariantLabelMapper = true;
    }

}