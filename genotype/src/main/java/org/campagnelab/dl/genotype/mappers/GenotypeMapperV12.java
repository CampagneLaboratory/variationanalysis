package org.campagnelab.dl.genotype.mappers;

/**
 * combined layer+inverse.
 */
public class GenotypeMapperV12 extends GenotypeMapperV11 {


    public GenotypeMapperV12() {
        super();
        sortCounts = true;
        withDistinctAlleleCounts = false;
        withCombinedLayer = true;
        hasIsVariantLabelMapper = false;
    }

}