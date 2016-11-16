package org.campagnelab.dl.somatic.mappers;

/**
 * Same as V13, but with FractionDifferences4 instead of 3., to avoid division by zero when sum to counts is 0
 * in one sample.
 */
public class FeatureMapperV18 extends NamingConcatFeatureMapper {
    public FeatureMapperV18() {
        super(new SimpleFeatureCalculator(true), new IndelFeatures(),
                new ReadIndexFeaturesFix(), new FractionDifferences4(), new MagnitudeFeatures2()
        );
    }
}
