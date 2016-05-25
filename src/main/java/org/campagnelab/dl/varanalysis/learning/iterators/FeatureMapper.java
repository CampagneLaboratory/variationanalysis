package org.campagnelab.dl.varanalysis.learning.iterators;

import org.campagnelab.dl.varanalysis.protobuf.BaseInformationRecords;
import org.nd4j.linalg.api.ndarray.INDArray;

/**
 * Created by fac2003 on 5/24/16.
 */
public interface FeatureMapper {

    /**
     * Return the number of features that this calculator will produce for each record.
     * The number of features is a constant across the datasets returned by the calculator.
     *
     * @return The number of features.
     */
    int numberOfFeatures();


    /**
     * Fill in features into the dataset, as index
     *  @param record        The record to convert to features & labels.
     * @param inputs        The features
     * @param indexOfRecord Index of the record in the destination dataset.
     */
    void mapFeatures(BaseInformationRecords.BaseInformationOrBuilder record, INDArray inputs, int indexOfRecord);

    /**
     * Produce the value of a given feature for the specified record.
     *
     * @param record       The record of interest.
     * @param featureIndex The index of the feature to produce/calculate
     * @return The value of the feature.
     */
    float produceFeature(BaseInformationRecords.BaseInformationOrBuilder record, int featureIndex);
}
