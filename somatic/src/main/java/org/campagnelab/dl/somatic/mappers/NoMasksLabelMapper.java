package org.campagnelab.dl.somatic.mappers;

import org.campagnelab.dl.framework.mappers.LabelMapper;
import org.campagnelab.dl.framework.mappers.MappedDimensions;
import org.nd4j.linalg.api.ndarray.INDArray;

/**
 * LabelMapper with no masks.
 * Created by fac2003 on 11/12/16.
 */
public abstract class NoMasksLabelMapper<RecordType> implements LabelMapper<RecordType> {
    /**
     * The default implementation returns a 1 dimension.
     * @return
     */
    @Override
    public MappedDimensions dimensions() {
        return new MappedDimensions(numberOfLabels());
    }

    @Override
    public boolean hasMask() {
        return false;
    }

    @Override
    public void maskLabels(RecordType record, INDArray mask, int indexOfRecord) {

    }

    @Override
    public boolean isMasked(RecordType record, int featureIndex) {
        return false;
    }

    @Override
    public void prepareToNormalize(RecordType record, int indexOfRecord) {

    }
    protected int sampleIndex;

    /**
     * Set the sample index on this mapper before calling configure.
     * @param sampleIndex index of the sample in the record whose features will be mapped.
     */
    public void setSampleIndex(int sampleIndex) {
        this.sampleIndex = sampleIndex;
    }
}
