package org.campagnelab.dl.genotype.mappers;

/**
 * Combined layer+inverse+isVariant label
 */
public class GenotypeMapperV13 extends GenotypeMapperV11 {


    public GenotypeMapperV13() {
        super();
        sortCounts = true;
        withDistinctAlleleCounts = false;
        withCombinedLayer = true;
        hasIsVariantLabelMapper = true;
    }

}