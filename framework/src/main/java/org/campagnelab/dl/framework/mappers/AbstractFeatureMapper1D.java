package org.campagnelab.dl.framework.mappers;

import org.nd4j.linalg.api.ndarray.INDArray;

/**
 * Stripped down version of AbstractFeatureMapper that implements some basic functions, should only be extended by intermediary mappers
 *
 * @author Remi Torracinta
 */
public abstract class AbstractFeatureMapper1D<RecordType> implements FeatureNameMapper<RecordType> {
    private int[] indices = new int[]{0,0};

    /**
     * This is a slow implementation for intermerdiary mappers to fulfill the requirements of the interface,
     * expecting that produceFeature will be used directly with these mappers.
     * Shoul not be used for a final mappers, which should override with a faster implementation.
     *
     * @param record        The record to convert to features & labels.
     * @param inputs        The features
     * @param indexOfRecord Index of the record in the destination dataset.
     */
    public void mapFeatures(RecordType record, INDArray inputs, int indexOfRecord) {


        for (int featureIndex = 0; featureIndex < numberOfFeatures(); featureIndex++) {
            indices[0] = indexOfRecord;
            indices[1] = featureIndex;
            inputs.putScalar(indices, produceFeature(record, featureIndex));
        }
    }

    public boolean hasMask() {
        return false;
    }

    public void maskFeatures(RecordType record, INDArray mask, int indexOfRecord) {
        // do nothing implementation.
    }

    public boolean isMasked(RecordType record, int featureIndex) {
        return false;
    }

    /**
     * Default implementation assumes a 1-d tensor.
     *
     * @return
     */
    @Override
    public MappedDimensions dimensions() {
        return new MappedDimensions(numberOfFeatures());
    }

}
