package org.campagnelab.dl.somatic.mappers;

/**
 * The FeatureMapper to test for the fourth iteration.
 * Created by rct66 on 5/31/16.
 */
public class FeatureMapperFT extends ConcatFeatureMapper {
    public FeatureMapperFT() {
        super(new FractionDifferences()
        );
    }
}
