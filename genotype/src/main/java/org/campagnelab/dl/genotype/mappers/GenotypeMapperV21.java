package org.campagnelab.dl.genotype.mappers;

/**
 * V19 (distinct alleles) with hasIsVariantLabelMapper.
 */
public class GenotypeMapperV21 extends GenotypeMapperV19 {


    public GenotypeMapperV21() {
        super();
        sortCounts = true;
        withDistinctAlleleCounts = true;
        withCombinedLayer = false;
        hasIsVariantLabelMapper = true;
    }

}