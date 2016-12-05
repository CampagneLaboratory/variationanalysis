package org.campagnelab.dl.somatic.mappers;

import org.campagnelab.dl.framework.iterators.ConcatFeatureMapper;

/**
 * Another feature mapper.
 *
 */
public class FeatureMapperV8 extends ConcatFeatureMapper {
    public FeatureMapperV8() {
        super(new SimpleFeatureCalculator(true),
                new SortedGenotypeAgreementMapper(),
                //new MagnitudeFeatures(),
               // new QualityFeatures(),
                new ReadIndexFeatures(),
                new FractionDifferences2()

        );
    }
}
