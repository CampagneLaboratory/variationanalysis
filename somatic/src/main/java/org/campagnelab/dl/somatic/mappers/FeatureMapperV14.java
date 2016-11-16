package org.campagnelab.dl.somatic.mappers;

/**
 * The FeatureMapper to test for the fourth iteration.
 * Created by rct66 on 5/31/16.
 */
public class FeatureMapperV14 extends ConcatFeatureMapper {
    public FeatureMapperV14() {
        super(new SimpleFeatureCalculator(true), new IndelFeatures(),
                new ReadIndexFeatures(), new FractionDifferences3(), new MagnitudeFeatures3()
        );
    }
}
