package org.campagnelab.dl.model.utils.mappers;

/**
 * Same as V13, but with FractionDifferences4 instead of 3., to avoid division by zero when sum to counts is 0
 * in one sample.
 */
public class FeatureMapperV16 extends ConcatFeatureMapper {
    public FeatureMapperV16() {
        super(new SimpleFeatureCalculator(true), new IndelFeatures(),
                new ReadIndexFeatures(), new FractionDifferences4(), new MagnitudeFeatures2()
        );
    }
}
