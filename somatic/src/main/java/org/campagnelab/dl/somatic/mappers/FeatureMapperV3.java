package org.campagnelab.dl.somatic.mappers;

/**
 * The FeatureMapper to test for the second iteration.
 * Created by rct66 on 5/31/16.
 */
public class FeatureMapperV3 extends ConcatFeatureMapper {
    public FeatureMapperV3() {
        super(new SimpleFeatureCalculator(), new MagnitudeFeatures(), new QualityFeatures());
    }
}
