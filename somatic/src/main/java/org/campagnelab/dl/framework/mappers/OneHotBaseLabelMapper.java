package org.campagnelab.dl.framework.mappers;

import org.nd4j.linalg.api.ndarray.INDArray;

import java.util.function.Function;

/**
 * Created by joshuacohen on 11/21/16.
 */
public class OneHotBaseLabelMapper<RecordType> implements LabelMapper<RecordType> {
    private int baseIndex;
    private int numLabels;
    private Function<RecordType, int[]> recordToLabel;

    private static final int[] indices = new int[]{0, 0};

    public OneHotBaseLabelMapper(int baseIndex, int numLabels, Function<RecordType, int[]> recordToLabel) {
        this.baseIndex = baseIndex;
        this.numLabels = numLabels;
        this.recordToLabel = recordToLabel;
    }

    @Override
    public void prepareToNormalize(RecordType record, int indexOfRecord) {

    }

    @Override
    public int numberOfLabels() {
        return numLabels;
    }

    @Override
    public MappedDimensions dimensions() {
        return new MappedDimensions(numberOfLabels());
    }

    @Override
    public void mapLabels(RecordType record, INDArray labels, int indexOfRecord) {
        indices[0] = indexOfRecord;
        for (int featureIndex = 0; featureIndex < numberOfLabels(); featureIndex++) {
            indices[1] = featureIndex;
            labels.putScalar(indices, produceLabel(record, featureIndex));
        }
    }

    @Override
    public float produceLabel(RecordType record, int labelIndex) {
        return labelIndex == recordToLabel.apply(record)[baseIndex] ? 1F : 0F;
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
}
