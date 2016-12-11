package org.campagnelab.dl.somatic.mappers;

import org.campagnelab.dl.framework.mappers.ConcatFeatureMapper;

/**
 * Same as V16, but with GenomicPositionMapper
 */
public class FeatureMapperV17 extends ConcatFeatureMapper {
    public FeatureMapperV17() {
        super(new SimpleFeatureCalculator(true), new IndelFeatures(),
                new ReadIndexFeatures(), new FractionDifferences4(), new MagnitudeFeatures2(),
                new GenomicPositionMapper()
        );
    }
}
