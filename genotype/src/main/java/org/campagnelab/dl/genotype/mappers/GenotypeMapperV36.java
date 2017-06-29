package org.campagnelab.dl.genotype.mappers;

/**
 * V35 (combined) with hasIsVariantLabelMapper. (Upgrade from V20FixRef).
 */
public class GenotypeMapperV36 extends GenotypeMapperV35 {

    public GenotypeMapperV36() {
        super(0);
    }

    public GenotypeMapperV36(int sampleIndex) {
        super(sampleIndex);
        sortCounts = true;
        withDistinctAlleleCounts = false;
        withCombinedLayer = true;
        withCombinedLayerRef = true;
        hasIsVariantLabelMapper = true;
    }
}

