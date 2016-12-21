package org.campagnelab.dl.framework.mappers;

/**
 * This mapper normalizes feature value by the absolute value of the maximum value obtained across all the delegate's
 * features. It accepts a single delegate, and normalizes it by dividing all the features by magnitude of the biggest feature.
 *
 * @author Remi Torracinta
 */

public class MaxNormalizationMapper<RecordType> extends AbstractFeatureMapper1D<RecordType> {

    FeatureNameMapper delegate;
    float absMax=Float.NEGATIVE_INFINITY;

    public MaxNormalizationMapper(FeatureNameMapper delegate) {
        this.delegate = delegate;
    }

    @Override
    public int numberOfFeatures() {
        return delegate.numberOfFeatures();
    }

    public void prepareToNormalize(RecordType record, int indexOfRecord) {
        for (int i = 0; i < numberOfFeatures(); i++) {
            absMax = Math.max(Math.abs(delegate.produceFeature(record, i)), absMax);
        }
    }

    public float produceFeature(RecordType record, int featureIndex) {
        return normalize(produceFeatureInternal(record, featureIndex), absMax);
    }

    @Override
    public String getFeatureName(int featureIndex) {
        return delegate.getFeatureName(featureIndex);
    }


    private float normalize(float value, float normalizationFactor) {
        if (normalizationFactor == 0) {
            return 0;
        }
        float normalized = value / normalizationFactor;
        assert normalized >= -1 && normalized <= 1 : "value must be normalized: " + normalized;
        return normalized;
    }


    private float produceFeatureInternal(RecordType record, int featureIndex) {
        return delegate.produceFeature(record, featureIndex);
    }

}


