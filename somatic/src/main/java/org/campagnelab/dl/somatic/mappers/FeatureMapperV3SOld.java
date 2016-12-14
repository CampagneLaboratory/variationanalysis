package org.campagnelab.dl.somatic.mappers;

import org.campagnelab.dl.framework.mappers.ConcatFeatureMapper;

/**
 * The FeatureMapper to test for the second iteration.
 * Created by rct66 on 5/31/16.
 */
public class FeatureMapperV3SOld extends ConcatFeatureMapper {
    public FeatureMapperV3SOld() {
        super( new SimpleFeatureCalculator(true), new SimpleFeatureCalculator(false) );
    }
}
