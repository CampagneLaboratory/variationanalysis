package org.campagnelab.dl.framework.mappers;

import org.campagnelab.dl.varanalysis.protobuf.BaseInformationRecords;
import org.nd4j.linalg.api.ndarray.INDArray;

/**
 * A wrapper to provide a name to a feature.
 * Created by fac2003 on 2/22/17.
 */
public abstract class NamedWrapper<RecordType> implements FeatureNameMapper<RecordType>, FeatureMapper<RecordType> {
    public NamedWrapper(FeatureMapper<RecordType> delegate) {
        this.delegate=delegate;
    }

    @Override
    public int numberOfFeatures() {
        return delegate.numberOfFeatures();
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
        delegate.mapFeatures(record, inputs, indexOfRecord);
    }

    @Override
    public boolean hasMask() {
        return delegate.hasMask();
    }

    @Override
    public void maskFeatures(RecordType record, INDArray mask, int indexOfRecord) {
        delegate.maskFeatures(record, mask, indexOfRecord);
    }

    @Override
    public boolean isMasked(RecordType record, int featureIndex) {
        return delegate.isMasked(record, featureIndex);
    }

    @Override
    public float produceFeature(RecordType record, int featureIndex) {
        return delegate.produceFeature(record, featureIndex);
    }

    FeatureMapper<RecordType> delegate;
}
