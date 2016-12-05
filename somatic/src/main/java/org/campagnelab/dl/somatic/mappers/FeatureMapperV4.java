package org.campagnelab.dl.somatic.mappers;

import org.campagnelab.dl.framework.iterators.ConcatFeatureMapper;

/**
 * The FeatureMapper to test for the fourth iteration.
 * Created by rct66 on 5/31/16.
 */
public class FeatureMapperV4 extends ConcatFeatureMapper {
    public FeatureMapperV4() {
        super(new SimpleFeatureCalculator(true),
                new SortedGenotypeAgreementMapper(),
                new MagnitudeFeatures(),
                new QualityFeatures(),
                new ReadIndexFeatures()

        );
    }
}
