package org.campagnelab.dl.genotype.mappers;

/**
 * V17 (combined) with hasIsVariantLabelMapper.
 */
public class GenotypeMapperV20FixRef extends GenotypeMapperV19 {


    public GenotypeMapperV20FixRef() {
        super();
        sortCounts = true;
        withDistinctAlleleCounts = false;
        withCombinedLayer = true;
        withCombinedLayerRef = true;
        hasIsVariantLabelMapper = true;
    }

}