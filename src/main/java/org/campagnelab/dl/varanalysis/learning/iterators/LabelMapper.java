package org.campagnelab.dl.varanalysis.learning.iterators;

import org.campagnelab.dl.varanalysis.protobuf.BaseInformationRecords;
import org.nd4j.linalg.api.ndarray.INDArray;

/**
 * Created by fac2003 on 5/24/16.
 */
public interface LabelMapper {
    /**
     * Return the number of labels that this calculator will produce for each record.
     * The number of labels is a constant across the datasets returned by the calculator.
     *
     * @return The number of labels.
     */
    int numberOfLabels();
    /**
     * Fill in features into the dataset, as index
     *  @param record        The record to convert to features & labels.

     * @param labels        The labels
     * @param indexOfRecord Index of the record in the destination dataset.
     */
    void mapLabels(BaseInformationRecords.BaseInformationOrBuilder record, INDArray labels, int indexOfRecord);
    /**
     * Produce the value of a given label for the specified record.
     *
     * @param record     The record of interest.
     * @param labelIndex The index of the label to produce/calculate
     * @return The value of the label.
     */
    float produceLabel(BaseInformationRecords.BaseInformationOrBuilder record, int labelIndex);
}
