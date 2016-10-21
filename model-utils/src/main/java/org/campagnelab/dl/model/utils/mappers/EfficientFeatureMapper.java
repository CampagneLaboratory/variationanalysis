package org.campagnelab.dl.model.utils.mappers;

import org.campagnelab.dl.varanalysis.protobuf.BaseInformationRecords;
import org.campagnelab.goby.baseinfo.SequenceBaseInformationReader;
import org.nd4j.linalg.api.ndarray.INDArray;

/**
 * A more efficient feature mapper interface, which does not depend on DL4J's INDArrays and the putScalar method.
 * Also adds a configure method to customize mapping for properties of entire datasets.
 *
 * @author Fabien Campagne Created by fac2003 on 10/19/16.
 *         New feature mappers are written with this interface, which must be used in place of FeatureMapper. When the interface
 *         is used, the more efficient data loading approach is used. Otherwise we fall back to the previous loading approach.
 * @since 1.0.3
 */
public interface EfficientFeatureMapper {
    /**
     * Configure the feature mapper for a specific set of sbi files. This method may access the properties of the reader
     * to retrieve statistics about the data being mapped to features (such as dataset wide normalization data).
     * @param reader
     */
    void configure(SequenceBaseInformationReader reader);
    /**
     * Return the number of features that this calculator will produce for each record.
     * The number of features is a constant across the datasets returned by the calculator.
     *
     * @return The number of features.
     */
    int numberOfFeatures();

    /**
     * This method must be called before produceFeature, so that the mapper has a chance to estimate the normalization factor.
     *
     * @param record        The record to convert to features & labels.
     * @param indexOfRecord Index of the record in the destination dataset.
     */
    public void prepareToNormalize(BaseInformationRecords.BaseInformationOrBuilder record, int indexOfRecord);

    /**
     * Produce the value of a given feature for the specified record. Will return a normalized feature
     * if the mapper implements normalization.
     *
     * @param record       The record of interest.
     * @param featureIndex The index of the feature to produce/calculate
     * @return The value of the feature.
     */
    float produceFeature(BaseInformationRecords.BaseInformationOrBuilder record, int featureIndex);

    /**
     * Fill in features into the dataset, as index
     *
     * @param record        The record to convert to features & labels.
     * @param inputs        The features will be written to this array. Must have size numberOfFeatures().
     * @param offset        The offset into the inputs array where features will be written consecutively.
     * @param indexOfRecord Index of the record in the destination dataset.
     */
    void mapFeatures(BaseInformationRecords.BaseInformationOrBuilder record, float[] inputs, int offset, int indexOfRecord);

}
