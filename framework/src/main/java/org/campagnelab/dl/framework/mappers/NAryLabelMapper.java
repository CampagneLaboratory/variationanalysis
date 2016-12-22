package org.campagnelab.dl.framework.mappers;

import org.nd4j.linalg.api.ndarray.INDArray;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

/**
 * Maps a one dimensional label of numLabels ones or zeros, with configuration over
 * whether the label should be masked or not
 *
 * Created by joshuacohen on 12/14/16.
 */
public class NAryLabelMapper<RecordType> implements LabelMapper<RecordType> {
    static private Logger LOG = LoggerFactory.getLogger(NAryLabelMapper.class);

    private int numLabels;
    private boolean isZero;
    private boolean hasMask;
    private boolean isMasked;
    private MappedDimensions dim;

    private int[] mapperIndices = new int[]{0, 0};

    /**
     * @param numLabels Number of individual labels the label mapper maps
     * @param isZero If true, label filled with zeros; if false, filled with ones
     * @param hasMask Whether or not the mapper has a mask
     * @param isMasked If hasMask is true, specifies if mask should be 1 (i.e., isMasked true) or 0
     */
    public NAryLabelMapper(int numLabels, boolean isZero, boolean hasMask, boolean isMasked) {
        this.numLabels = numLabels;
        this.isZero = isZero;
        if (numLabels > 1 && !isZero) {
            LOG.warn("Mapping multiple 1s may break one-hot encoding");
        }
        this.hasMask = hasMask;
        this.isMasked = isMasked;
        dim = new MappedDimensions(numLabels);
    }

    @Override
    public int numberOfLabels() {
        return numLabels;
    }

    @Override
    public MappedDimensions dimensions() {
        return dim;
    }

    @Override
    public void prepareToNormalize(RecordType record, int indexOfRecord) {

    }

    @Override
    public void mapLabels(RecordType record, INDArray inputs, int indexOfRecord) {
        mapperIndices[0] = indexOfRecord;
        for (int i = 0; i < numLabels; i++) {
            mapperIndices[1] = i;
            inputs.putScalar(mapperIndices, produceLabel(record, i));
        }
    }

    @Override
    public boolean hasMask() {
        return hasMask;
    }

    @Override
    public void maskLabels(RecordType record, INDArray mask, int indexOfRecord) {
        mask.putScalar(indexOfRecord, isMasked(record, -1) ? 1F : 0F);
    }

    @Override
    public boolean isMasked(RecordType record, int labelIndex) {
        return isMasked;
    }

    @Override
    public float produceLabel(RecordType record, int labelIndex) {
        return isZero ? 0F : 1F;
    }
}

