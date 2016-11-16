package org.campagnelab.dl.framework.mappers;

import org.nd4j.linalg.api.ndarray.INDArray;

/**
 * Created by fac2003 on 5/24/16.
 */
public interface LabelMapper<RecordType> {
    /**
     * Return the number of labels that this calculator will produce for each record.
     * The number of labels is a constant across the datasets returned by the calculator.
     *
     * @return The number of labels.
     */
    int numberOfLabels();

    /**
     * Fill in features into the dataset, as index
     *
     * @param record        The record to convert to features & labels.
     * @param labels        The labels
     * @param indexOfRecord Index of the record in the destination dataset.
     */
    void mapLabels(RecordType record, INDArray labels, int indexOfRecord);

    /**
     * Produce the value of a given label for the specified record.
     *
     * @param record     The record of interest.
     * @param labelIndex The index of the label to produce/calculate
     * @return The value of the label.
     */
    float produceLabel(RecordType record, int labelIndex);

    /**
     * Return true if the mapper creates an input mask (maskFeatures is implemented).
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
}
