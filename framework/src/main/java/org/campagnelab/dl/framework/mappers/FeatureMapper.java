package org.campagnelab.dl.framework.mappers;

import org.nd4j.linalg.api.ndarray.INDArray;

/**
 * FeatureMapper instances convert records to mapped features suitable to train a neural net or computation graph.
 * Created by fac2003 on 5/24/16.
 */
public interface FeatureMapper<RecordType> {

    /**
     * Return the number of features that this calculator will produce for each record.
     * The number of features is a constant across the datasets returned by the calculator.
     *
     * @return The number of features.
     */
    int numberOfFeatures();

    /**
     * Dimensions of the feature and mask tensors. These dimensions must match the input INDArray used to call
     * mapFeatures and maskLabels.
     * @return mapped tensor dimensions.
     */
    MappedDimensions dimensions();

    /**
     * This method must be called before produceFeature, so that the mapper has a chance to estimate the normalization factor.
     *
     * @param record        The record to convert to features & labels.
     * @param indexOfRecord Index of the record in the destination dataset.
     */
    public void prepareToNormalize(RecordType record, int indexOfRecord);

    /**
     * Fill in features into the dataset, as index
     *
     * @param record        The record to convert to features & labels.
     * @param inputs        The features
     * @param indexOfRecord Index of the record in the destination dataset.
     */
    void mapFeatures(RecordType record, INDArray inputs, int indexOfRecord);

    /**
     * Return true if the mapper creates an input mask (maskLabels is implemented).
     *
     * @return
     */
    boolean hasMask();

    /**
     * Fill in the feature mask. The method is only called if hasMask returns true.
     *
     * @param record        The record to convert to features & labels.
     * @param mask          The feature mask
     * @param indexOfRecord Index of the record in the destination dataset.
     */
    void maskFeatures(RecordType record, INDArray mask, int indexOfRecord);

    /**
     * Determine if a feature needs to be masked or not. The method is only called if hasMask returns true.
     *
     * @param record       the record for which masking may be needed.
     * @param featureIndex index of the features that may need masking
     * @return True if feature must be masked, false otherwise.
     */
    boolean isMasked(RecordType record, int featureIndex);

    /**
     * Produce the value of a given feature for the specified record. Will return a normalized feature
     * if the mapper implements normalization.
     *
     * @param record       The record of interest.
     * @param featureIndex The index of the feature to produce/calculate
     * @return The value of the feature.
     */
    float produceFeature(RecordType record, int featureIndex);


}
