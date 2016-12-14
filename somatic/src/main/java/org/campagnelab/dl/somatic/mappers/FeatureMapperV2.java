package org.campagnelab.dl.somatic.mappers;

import org.campagnelab.dl.framework.mappers.ConcatFeatureMapper;

/**
 * The FeatureMapper to test for the second iteration.
 * Created by fac2003 on 5/24/16.
 */
public class FeatureMapperV2 extends ConcatFeatureMapper {
    public FeatureMapperV2() {
        super(new SimpleFeatureCalculator(true), new MagnitudeFeatures());
    }
}
