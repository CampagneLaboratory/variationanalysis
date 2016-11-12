package org.campagnelab.dl.model.utils.mappers;

import org.nd4j.linalg.api.ndarray.INDArray;

/**
 * LabelMapper with no masks.
 * Created by fac2003 on 11/12/16.
 */
public abstract class NoMasksLabelMapper<RecordType> implements LabelMapper<RecordType> {
    @Override
    public boolean hasMask() {
        return false;
    }

    @Override
    public void maskFeatures(RecordType record, INDArray mask, int indexOfRecord) {

    }

    @Override
    public boolean isMasked(RecordType record, int featureIndex) {
        return false;
    }
}
