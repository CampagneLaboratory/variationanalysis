package org.campagnelab.dl.varanalysis.mappers;

import org.nd4j.linalg.api.ndarray.INDArray;

/**
 * Created by fac2003 on 11/11/16.
 */
public abstract class NoMaskFeatureMapper<RecordType> implements FeatureMapper<RecordType> {
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
