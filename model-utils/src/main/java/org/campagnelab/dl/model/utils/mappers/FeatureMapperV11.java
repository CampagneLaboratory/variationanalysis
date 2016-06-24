package org.campagnelab.dl.model.utils.mappers;

/**
 * The FeatureMapper to test for the fourth iteration.
 * Created by rct66 on 5/31/16.
 */
public class FeatureMapperV11 extends ConcatFeatureMapper {
    public FeatureMapperV11() {
        super(new SimpleFeatureCalculator(false), new SimpleFeatureCalculator(true), new MagnitudeFeatures(), new QualityFeatures(),
                new ReadIndexFeatures(), new FractionDifferences2()
        );
    }
}
