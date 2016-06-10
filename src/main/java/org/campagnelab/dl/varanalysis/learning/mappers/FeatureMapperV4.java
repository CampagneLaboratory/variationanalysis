package org.campagnelab.dl.varanalysis.learning.mappers;

/**
 * The FeatureMapper to test for the fourth iteration.
 * Created by rct66 on 5/31/16.
 */
public class FeatureMapperV4 extends ConcatFeatureMapper {
    public FeatureMapperV4() {
        super(new MagnitudeFeatures(), new SimpleFeatureCalculator(false),
                new MagnitudeFeatures(), new SimpleFeatureCalculator(true),new QualityFeatures(),
                new ReadIndexFeatures()
        );
    }
}
