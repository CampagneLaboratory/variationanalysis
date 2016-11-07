package org.campagnelab.dl.model.utils.mappers;

import org.campagnelab.dl.varanalysis.protobuf.BaseInformationRecords;
import org.nd4j.linalg.api.ndarray.INDArray;

/**
 * A more efficient label mapper interface.
 */
public interface EfficientLabelMapper<RecordType> {
    /**
     * Return the number of labels that this calculator will produce for each record.
     * The number of labels is a constant across the datasets returned by the calculator.
     *
     * @return The number of labels.
     */
    int numberOfLabels();

    /**
     * Produce labels for a dataset.
     *
     * @param record        The record to convert to features & labels.
     * @param labels        The labels will be written to this array. Must have size numberOfLabels().
     * @param offset        into the labels array.
     * @param indexOfRecord Index of the record in the destination dataset.
     */
    void mapLabels(RecordType record, float[] labels, int offset, int indexOfRecord);

    /**
     * Produce the value of a given label for the specified record.
     *
     * @param record     The record of interest.
     * @param labelIndex The index of the label to produce/calculate
     * @return The value of the label.
     */
    float produceLabel(RecordType record, int labelIndex);
}
