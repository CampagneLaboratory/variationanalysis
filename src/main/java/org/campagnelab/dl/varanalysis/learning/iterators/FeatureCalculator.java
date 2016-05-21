package org.campagnelab.dl.varanalysis.learning.iterators;

import org.campagnelab.dl.varanalysis.format.PosRecord;
import org.nd4j.linalg.api.ndarray.INDArray;
import org.nd4j.linalg.dataset.DataSet;

/**
 * Maps base information to features for maching learning with neural nets.
 * Created by fac2003 on 5/21/16.
 *
 * @author Fabien Campagne
 */
interface FeatureCalculator {

    /**
     * Return the number of features that this calculator will produce for each record.
     * The number of features is a constant across the datasets returned by the calculator.
     *
     * @return The number of features.
     */
    int numberOfFeatures();

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
     * @param inputs        The features
     * @param labels        The labels
     * @param indexOfRecord Index of the record in the destination dataset.
     */
    void map(PosRecord record, INDArray inputs, INDArray labels, int indexOfRecord);

    /**
     * Produce the value of a given feature for the specified record.
     *
     * @param record       The record of interest.
     * @param featureIndex The index of the feature to produce/calculate
     * @return The value of the feature.
     */
    float produceFeature(PosRecord record, int featureIndex);

    /**
     * Produce the value of a given label for the specified record.
     *
     * @param record     The record of interest.
     * @param labelIndex The index of the label to produce/calculate
     * @return The value of the label.
     */
    float produceLabel(PosRecord record, int labelIndex);
}
