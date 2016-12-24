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