package org.campagnelab.dl.model.utils.mappers;

/**
 * The FeatureMapper to test for the fourth iteration.
 * Created by rct66 on 5/31/16.
 */
public class FeatureMapperV13 extends ConcatFeatureMapper {
    public FeatureMapperV13() {
        super(new SimpleFeatureCalculator(true), new IndelFeatures(),
                new ReadIndexFeatures(), new FractionDifferences3(), new MagnitudeFeatures2()
        );
    }
}
