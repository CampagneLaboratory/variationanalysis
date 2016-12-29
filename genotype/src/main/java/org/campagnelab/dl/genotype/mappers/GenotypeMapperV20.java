package org.campagnelab.dl.genotype.mappers;

/**
 * V17 (combined) with hasIsVariantLabelMapper.
 */
public class GenotypeMapperV20 extends GenotypeMapperV19 {


    public GenotypeMapperV20() {
        super();
        sortCounts = true;
        withDistinctAlleleCounts = false;
        withCombinedLayer = true;
        hasIsVariantLabelMapper = true;
    }

}