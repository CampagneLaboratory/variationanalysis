package org.campagnelab.dl.varanalysis.learning.mappers;

/**
 * The FeatureMapper to test for the fourth iteration.
 * Created by rct66 on 5/31/16.
 */
public class FeatureMapperVT extends ConcatFeatureMapper {
    public FeatureMapperVT() {
        super(new SimpleFeatureCalculator(false), new SimpleFeatureCalculator(true), new MagnitudeFeatures(), new QualityFeatures(),
                new ReadIndexFeatures()
        );
    }
}
