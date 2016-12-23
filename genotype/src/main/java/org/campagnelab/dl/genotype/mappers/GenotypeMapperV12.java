package org.campagnelab.dl.genotype.mappers;

import org.campagnelab.dl.framework.mappers.FeatureNameMapper;
import org.campagnelab.dl.somatic.mappers.DensityMapper;
import org.campagnelab.dl.somatic.mappers.GenomicContextMapper;
import org.campagnelab.dl.somatic.mappers.NamingConcatFeatureMapper;
import org.campagnelab.dl.somatic.mappers.functional.TraversalHelper;
import org.campagnelab.dl.varanalysis.protobuf.BaseInformationRecords;
import org.nd4j.linalg.api.ndarray.INDArray;

import java.util.Properties;

/**
 * Based of V10,  switch from Max normalization to inverse normalization.
 */
public class GenotypeMapperV12 extends GenotypeMapperV4 {


    public GenotypeMapperV12() {
        super();
        sortCounts = true;
        withDistinctAlleleCounts = true;
        hasIsVariantLabelMapper = true;
    }

}