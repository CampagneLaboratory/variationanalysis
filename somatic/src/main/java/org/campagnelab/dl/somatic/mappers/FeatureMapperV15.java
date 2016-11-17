package org.campagnelab.dl.somatic.mappers;

import org.campagnelab.dl.framework.iterators.ConcatFeatureMapper;

/**
 * The FeatureMapper to test for the fourth iteration.
 * Created by rct66 on 5/31/16.
 */
public class FeatureMapperV15 extends ConcatFeatureMapper {
    public FeatureMapperV15() {
        super(new SimpleFeatureCalculator(true), new IndelFeatures(),
                new ReadIndexFeatures(), new FractionDifferences3(), new MagnitudeFeatures2(),
                new GenomicPositionMapper()
        );
    }
}
