package org.campagnelab.dl.framework.mappers;

import org.nd4j.linalg.api.ndarray.INDArray;

/**
 * Feature mapper created from a delegate label mapper, with all of the relevant
 * feature methods (produceFeature, mapFeatures, numberOfFeatures) corresponding to the same
 * methods in label (i.e. produceLabel, mapLabels, numberOfLabels)
 *
 * Created by joshuacohen on 12/14/16.
 */
public class FeatureFromLabelMapper<RecordType> implements FeatureMapper<RecordType> {
    LabelMapper<RecordType> delegate;

    /**
     * Creates a feature mapper from a label mapper
     * @param delegate label mapper delegate
     */
    public FeatureFromLabelMapper(LabelMapper<RecordType> delegate) {
        this.delegate = delegate;
    }

    @Override
    public int numberOfFeatures() {
        return delegate.numberOfLabels();
    }

    @Override
    public MappedDimensions dimensions() {
        return delegate.dimensions();
    }

    @Override
    public void prepareToNormalize(RecordType record, int indexOfRecord) {
        delegate.prepareToNormalize(record, indexOfRecord);
    }

    @Override
    public void mapFeatures(RecordType record, INDArray inputs, int indexOfRecord) {
        delegate.mapLabels(record, inputs, indexOfRecord);
    }

    @Override
    public boolean hasMask() {
        return delegate.hasMask();
    }

    @Override
    public void maskFeatures(RecordType record, INDArray mask, int indexOfRecord) {
        delegate.maskLabels(record, mask, indexOfRecord);
    }

    @Override
    public boolean isMasked(RecordType record, int featureIndex) {
        return delegate.isMasked(record, featureIndex);
    }

    @Override
    public float produceFeature(RecordType record, int featureIndex) {
        return delegate.produceLabel(record, featureIndex);
    }
}
