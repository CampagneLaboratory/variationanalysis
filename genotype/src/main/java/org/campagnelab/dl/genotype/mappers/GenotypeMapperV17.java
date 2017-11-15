package org.campagnelab.dl.genotype.mappers;

import org.campagnelab.dl.framework.mappers.FeatureNameMapper;
import org.campagnelab.dl.varanalysis.protobuf.BaseInformationRecords;

/**
 * V16 (=v15+distanceReads), + new PairFlag mappers
 * distinct alleles+inverse+ reduced number of bins for quality scores, mapping qual. Adding
 * bins for numVariationsInRead, targetAlignedLength, queryAlignedLength.
 */
public class GenotypeMapperV17 extends GenotypeMapperV16 {


    private FeatureNameMapper<BaseInformationRecords.BaseInformationOrBuilder> delegate;
    //default sampleIndex is zero, adjustable with setter
    private int sampleIndex = 0;

    public GenotypeMapperV17() {
        super();
        sortCounts = true;
        withDistinctAlleleCounts = true;
        withCombinedLayer = false;
        MAX_GENOTYPES = 3;
    }


}
