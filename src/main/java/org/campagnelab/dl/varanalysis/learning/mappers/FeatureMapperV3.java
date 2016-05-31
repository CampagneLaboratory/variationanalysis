package org.campagnelab.dl.varanalysis.learning.mappers;

import org.campagnelab.dl.varanalysis.learning.iterators.MagnitudeFeatures;
import org.campagnelab.dl.varanalysis.learning.iterators.QualityFeatures;
import org.campagnelab.dl.varanalysis.learning.iterators.SimpleFeatureCalculator;

/**
 * The FeatureMapper to test for the second iteration.
 * Created by rct66 on 5/31/16.
 */
public class FeatureMapperV3 extends ConcatFeatureMapper {
    public FeatureMapperV3() {
        super(new ConcatFeatureMapper(new SimpleFeatureCalculator(), new MagnitudeFeatures(), new QualityFeatures()));
    }
}
