package org.campagnelab.dl.somatic.mappers.trio;

import org.campagnelab.dl.somatic.mappers.NamingConcatFeatureMapper;

/**
 * Same as V13, but with FractionDifferences4 instead of 3., to avoid division by zero when sum to counts is 0
 * in one sample.
 */
public class FeatureMapperV18Trio extends NamingConcatFeatureMapper {
    public FeatureMapperV18Trio() {
        super(new SimpleFeatureCalculatorTrio(true), new IndelFeaturesTrio(),
                new ReadIndexFeaturesTrio(), new FractionDifferences4Trio(0), new FractionDifferences4Trio(1), new MagnitudeFeatures2Trio()
        );
    }
}
