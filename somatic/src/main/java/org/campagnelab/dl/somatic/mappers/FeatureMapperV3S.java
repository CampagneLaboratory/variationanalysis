package org.campagnelab.dl.somatic.mappers;

import org.campagnelab.dl.framework.mappers.ConcatFeatureMapper;

/**
 * The FeatureMapper to test for the second iteration.
 * Created by rct66 on 5/31/16.
 */
public class FeatureMapperV3S extends ConcatFeatureMapper {
    public FeatureMapperV3S() {
        super( new MagnitudeFeatures(), new SimpleFeatureCalculator(false),new MagnitudeFeatures(), new SimpleFeatureCalculator(true),new QualityFeatures() );
    }
}
