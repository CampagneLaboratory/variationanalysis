package org.campagnelab.dl.framework.mappers;

import it.unimi.dsi.fastutil.floats.FloatArrayList;
import org.nd4j.linalg.api.ndarray.INDArray;

/**
 * This is a normalizing mapper which divides feature values by their mean.
 *
 * @author Fabien Campagne
 */

public class MeanNormalizationMapper<RecordType> extends AbstractFeatureMapper1D<RecordType> {

    private final boolean dividebyStdev;
    FeatureNameMapper<RecordType> delegate;
    float mean = 0;
    boolean normalizedCalled = false;
    private double stdev;
    public MeanNormalizationMapper(FeatureNameMapper<RecordType> delegate) {
        this(delegate,false);
    }
    public MeanNormalizationMapper(FeatureNameMapper<RecordType> delegate, boolean dividebyStdev) {
        this.delegate = delegate;
        this.dividebyStdev = dividebyStdev;
    }

    @Override
    public int numberOfFeatures() {
        return delegate.numberOfFeatures();
    }


    @Override
    public void mapFeatures(RecordType record, INDArray inputs, int indexOfRecord) {
        super.mapFeatures(record, inputs, indexOfRecord);
        normalizedCalled = false;
    }

    private FloatArrayList values = new FloatArrayList();

    public void prepareToNormalize(RecordType record, int indexOfRecord) {
        mean = 0;
        stdev = 0;
        int count = 0;
        delegate.prepareToNormalize(record, indexOfRecord);
        values.clear();
        for (int i = 0; i < numberOfFeatures(); i++) {

            final float v = delegate.produceFeature(record, i);
            values.add(v);
            mean += v;
            count += 1;
        }
        mean /= count;
        double variance = 0;
        for (double value : values) {
            final double difference = value - mean;
            variance += difference * difference;
        }
        stdev = Math.sqrt(variance);

        delegate.prepareToNormalize(record, indexOfRecord);
        normalizedCalled = true;
    }


    public float produceFeature(RecordType record, int featureIndex) {
        assert normalizedCalled == true : "normalized must be called before produceFeature";
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
        float normalized = value - mean;
        if (dividebyStdev) {
            normalized /= stdev;

        }
        return normalized;
    }


    private float produceFeatureInternal(RecordType record, int featureIndex) {
        return delegate.produceFeature(record, featureIndex);
    }

}


