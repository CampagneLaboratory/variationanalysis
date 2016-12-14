package org.campagnelab.dl.somatic.mappers;

import org.campagnelab.dl.framework.mappers.ConcatFeatureMapper;

/**
 * The FeatureMapper to test for the fourth iteration.
 * Created by rct66 on 5/31/16.
 */
public class FeatureMapperV6 extends ConcatFeatureMapper {
    public FeatureMapperV6() {
        super(new FractionDifferences(), new MagnitudeFeatures(), new ReadIndexFeatures(), new QualityFeatures());
    }
}
