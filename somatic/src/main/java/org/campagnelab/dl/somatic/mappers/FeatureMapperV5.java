package org.campagnelab.dl.somatic.mappers;

import org.campagnelab.dl.framework.mappers.ConcatFeatureMapper;

/**
 * The FeatureMapper to test for the fourth iteration.
 * Created by rct66 on 5/31/16.
 */
public class FeatureMapperV5 extends ConcatFeatureMapper {
    public FeatureMapperV5() {
        super(new MagnitudeFeatures(), new SimpleFeatureCalculator(false), new QualityFeatures(), new ReadIndexFeatures());
    }
}
