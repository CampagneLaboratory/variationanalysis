package org.campagnelab.dl.genotype.mappers;

import org.campagnelab.dl.framework.mappers.FeatureNameMapper;
import org.campagnelab.dl.framework.mappers.MaxNormalizationMapper;
import org.campagnelab.dl.somatic.mappers.DensityMapper;
import org.campagnelab.dl.somatic.mappers.GenomicContextMapper;
import org.campagnelab.dl.somatic.mappers.NamingConcatFeatureMapper;
import org.campagnelab.dl.somatic.mappers.functional.TraversalHelper;
import org.campagnelab.dl.varanalysis.protobuf.BaseInformationRecords;
import org.nd4j.linalg.api.ndarray.INDArray;

import java.util.Properties;

/**
 * This mapper sorts counts.
 */
public class GenotypeMapperV5 extends GenotypeMapperV4 {
    public GenotypeMapperV5() {
        withDistinctAlleleCounts=true;
    }
}
