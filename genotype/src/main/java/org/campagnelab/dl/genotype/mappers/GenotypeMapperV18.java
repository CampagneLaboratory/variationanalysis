package org.campagnelab.dl.genotype.mappers;

/**
 * V17 (combined) with hasIsVariantLabelMapper.
 */
public class GenotypeMapperV18 extends GenotypeMapperV17 {


    public GenotypeMapperV18() {
        super();
        sortCounts = true;
        withDistinctAlleleCounts = false;
        withCombinedLayer = true;
        hasIsVariantLabelMapper = true;
    }

}