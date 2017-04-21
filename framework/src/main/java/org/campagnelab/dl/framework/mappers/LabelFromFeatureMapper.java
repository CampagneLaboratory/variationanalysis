package org.campagnelab.dl.framework.mappers;

import org.nd4j.linalg.api.ndarray.INDArray;

import java.util.Properties;

/**
 * Label mapper created from a delegate feature mapper, with all of the relevant
 * label methods (produceLabel, mapLabels, numberOfLabels) corresponding to the same
 * methods in feature (i.e. produceFeature, mapFeatures, numberOfFeatures)
 * <p>
 * Created by joshuacohen on 12/14/16.
 */
public class LabelFromFeatureMapper<RecordType> implements LabelMapper<RecordType>, ConfigurableLabelMapper {
    FeatureMapper<RecordType> delegate;

    /**
     * Creates a label mapper from a feature mapper
     *
     * @param delegate feature mapper delegate
     */
    public LabelFromFeatureMapper(FeatureMapper<RecordType> delegate) {
        this.delegate = delegate;
    }

    @Override
    public int numberOfLabels() {
        return delegate.numberOfFeatures();
    }

    @Override
    public MappedDimensions dimensions() {
        return delegate.dimensions();
    }

    @Override
    public void mapLabels(RecordType record, INDArray labels, int indexOfRecord) {
        delegate.mapFeatures(record, labels, indexOfRecord);
    }

    @Override
    public float produceLabel(RecordType record, int labelIndex) {
        return delegate.produceFeature(record, labelIndex);
    }

    @Override
    public boolean hasMask() {
        return delegate.hasMask();
    }

    @Override
    public void maskLabels(RecordType record, INDArray mask, int indexOfRecord) {
        delegate.maskFeatures(record, mask, indexOfRecord);
    }

    @Override
    public boolean isMasked(RecordType record, int featureIndex) {
        return delegate.isMasked(record, featureIndex);
    }

    @Override
    public void prepareToNormalize(RecordType record, int indexOfRecord) {
        delegate.prepareToNormalize(record, indexOfRecord);
    }

    @Override
    public void configure(Properties readerProperties) {
        if (delegate instanceof ConfigurableFeatureMapper) {
            ((ConfigurableFeatureMapper) delegate).configure(readerProperties);
        }
    }
}
