package org.campagnelab.dl.genotype.mappers;

/**
 * V22  withDistinctAlleleCounts  with hasIsVariantLabelMapper.
 */
public class GenotypeMapperV24 extends GenotypeMapperV22 {


    public GenotypeMapperV24() {
        super();
        sortCounts = true;
        withDistinctAlleleCounts = true;
        withCombinedLayer = false;
        hasIsVariantLabelMapper = true;
    }

}