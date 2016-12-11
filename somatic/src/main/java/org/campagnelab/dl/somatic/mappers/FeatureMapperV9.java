package org.campagnelab.dl.somatic.mappers;

import org.campagnelab.dl.framework.mappers.ConcatFeatureMapper;

/**
 * Another feature mapper.
 */
public class FeatureMapperV9 extends ConcatFeatureMapper {
    public FeatureMapperV9() {
        super(new SimpleFeatureCalculator(true),
                new SimpleFeatureCalculator(false),
                new SortedGenotypeAgreementMapper(),
                //new MagnitudeFeatures(),
                // new QualityFeatures(),
                new ReadIndexFeatures(),
                new FractionDifferences2()

        );
    }
}
