package org.campagnelab.dl.framework.mappers;

import org.nd4j.linalg.api.ndarray.INDArray;

/**
 * This is a normalizing mapper which divides feature values by their mean.
 * @author Fabien Campagne
 */

public class MeanNormalizationMapper<RecordType> extends AbstractFeatureMapper1D<RecordType> {

    FeatureNameMapper<RecordType> delegate;
    float mean=0;

    public MeanNormalizationMapper(FeatureNameMapper<RecordType> delegate) {
        this.delegate = delegate;
    }

    @Override
    public int numberOfFeatures() {
        return delegate.numberOfFeatures();
    }



    public void prepareToNormalize(RecordType record, int indexOfRecord) {
        mean=0;
        int count=0;
        delegate.prepareToNormalize(record, indexOfRecord);
        for (int i = 0; i < numberOfFeatures(); i++) {

            mean+=delegate.produceFeature(record, i);
            count+=1;
        }
        mean/=count;
    }


    public float produceFeature(RecordType record, int featureIndex) {
        return normalize(produceFeatureInternal(record, featureIndex), mean);
    }

    @Override
    public String getFeatureName(int featureIndex) {
        return delegate.getFeatureName(featureIndex);
    }


    private float normalize(float value, float normalizationFactor) {
        if (normalizationFactor == 0) {
            return 0;
        }
        return value / normalizationFactor;
        }


    private float produceFeatureInternal(RecordType record, int featureIndex) {
        return delegate.produceFeature(record, featureIndex);
    }

}


