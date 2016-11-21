package org.campagnelab.dl.genotype.mappers;

import it.unimi.dsi.fastutil.objects.ObjectArrayList;
import org.campagnelab.dl.framework.mappers.FeatureNameMapper;
import org.campagnelab.dl.somatic.genotypes.GenotypeCountFactory;
import org.campagnelab.dl.somatic.mappers.GenotypeCount;
import org.campagnelab.dl.varanalysis.protobuf.BaseInformationRecords;
import org.nd4j.linalg.api.ndarray.INDArray;

import java.util.Arrays;
import java.util.Collections;

/**
 * Stripped down version of AbstractFeatureMapper that implements some basic functions, should only be extended by intermediary mappers
 * @author Remi Torracinta
 */
public abstract class AbstractFeatureMapperStripped<T extends BaseInformationRecords.BaseInformationOrBuilder > implements FeatureNameMapper<T> {

    private float[] buffer;

    protected float[] getBuffer() {

        if (buffer==null) {
            buffer = new float[numberOfFeatures()];
        }else{
            Arrays.fill(buffer,0f);
        }
        return buffer;
    }



    int[] indices = new int[]{0, 0};


    /**
     * This is a slow implementation for intermerdiary mappers to fulfill the requirements of the interface,
     * expecting that produceFeature will be used directly with these mappers.
     * Shoul not be used for a final mappers, which should override with a faster implementation.
     * @param record        The record to convert to features & labels.
     * @param inputs        The features
     * @param indexOfRecord Index of the record in the destination dataset.
     */
    @Deprecated
    public void mapFeatures(T record, INDArray inputs, int indexOfRecord) {
        indices[0] = indexOfRecord;
        inputs.putScalar(indices,produceFeature(record,0));
    }


    /**
     * This is a slow implementation for intermerdiary mappers to fulfill the requirements of the interface,
     * expecting that produceFeature will be used directly with these mappers.
     * Shoul not be used for a final mappers, which should override with a faster implementation.
     * @param record
     * @param inputs
     * @param offset
     * @param indexOfRecord
     */
    @Deprecated
    public void mapFeatures(T record, float[] inputs, int offset, int indexOfRecord) {
        inputs[offset] = produceFeature(record,0);
    }


    public boolean hasMask() {
        return false;
    }

    public void maskFeatures(T record, INDArray mask, int indexOfRecord) {
        // do nothing implementation.
    }

    public boolean isMasked(T record, int featureIndex) {
        return false;
    }


}
