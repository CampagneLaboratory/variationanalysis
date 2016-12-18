package org.campagnelab.dl.framework.mappers;

/**
 * Created by rct2002 on 7/4/16.
 */
public interface FeatureNameMapper<T> extends FeatureMapper<T> {

    /**
     * Produce the name of a given feature.
     *
     * @param featureIndex The index of the feature to get the name of
     * @return The name of the feature.
     */
    String getFeatureName(int featureIndex);
}
