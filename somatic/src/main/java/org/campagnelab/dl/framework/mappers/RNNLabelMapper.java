package org.campagnelab.dl.framework.mappers;

import org.nd4j.linalg.api.ndarray.INDArray;

import java.util.function.Function;

/**
 * Created by joshuacohen on 11/21/16.
 */
public class RNNLabelMapper<RecordType> implements LabelMapper<RecordType> {
    private int labelsPerTimeStep;
    private OneHotBaseLabelMapper<RecordType>[] delegates;
    private Function<RecordType, Integer> recordToSequenceLength;

    private int[] indicesMapper = new int[]{0, 0, 0};
    private int[] indicesMasker = new int[]{0, 0};

    public RNNLabelMapper(int labelsPerTimeStep, int maxSequenceLength, Function<RecordType, int[]> recordToLabel,
                          Function<RecordType, Integer> recordToSequenceLength) {
        this.recordToSequenceLength = recordToSequenceLength;
        delegates = new OneHotBaseLabelMapper[maxSequenceLength];
        for (int i = 0; i < maxSequenceLength; i++) {
            delegates[i] = new OneHotBaseLabelMapper<>(i, labelsPerTimeStep, recordToLabel);
        }
        this.labelsPerTimeStep = delegates[0].numberOfLabels();
        for (OneHotBaseLabelMapper<RecordType> oneHotBaseMapper : delegates) {
            if (oneHotBaseMapper.numberOfLabels() != labelsPerTimeStep) {
                throw new RuntimeException("All delegate one hot base mappers should have same number of labels");
            }
        }
    }

    @Override
    public int numberOfLabels() {
        return delegates.length * labelsPerTimeStep;
    }

    @Override
    public MappedDimensions dimensions() {
        return new MappedDimensions(labelsPerTimeStep, delegates.length);
    }

    @Override
    public void mapLabels(RecordType record, INDArray labels, int indexOfRecord) {
        indicesMapper[0] = indexOfRecord;
        for (int i = 0; i < delegates.length; i++) {
            indicesMapper[2] = i;
            for (int j = 0; j < delegates[i].numberOfLabels(); j++) {
                indicesMapper[1] = j;
                if (i < recordToSequenceLength.apply(record)) {
                    labels.putScalar(indicesMapper, delegates[i].produceLabel(record, j));
                } else {
                    labels.putScalar(indicesMapper, 0F);
                }
            }
        }
    }

    @Override
    public boolean hasMask() {
        return true;
    }

    @Override
    public void maskLabels(RecordType record, INDArray mask, int indexOfRecord) {
        indicesMasker[0] = indexOfRecord;
        for (int i = 0; i < delegates.length; i++) {
            indicesMasker[1] = i;
            mask.putScalar(indicesMasker, i < recordToSequenceLength.apply(record) ? 1F : 0F);
        }
    }

    @Override
    public float produceLabel(RecordType record, int labelIndex) {
        int delegateIdx = labelIndex / labelsPerTimeStep;
        int labelInDelegateIdx = labelIndex % labelsPerTimeStep;
        return delegates[delegateIdx].produceLabel(record, labelInDelegateIdx);
    }

    @Override
    public boolean isMasked(RecordType record, int labelIndex) {
        int delegateIdx = labelIndex / labelsPerTimeStep;
        int labelInDelegateIdx = labelIndex % labelsPerTimeStep;
        return delegates[delegateIdx].isMasked(record, labelInDelegateIdx);
    }
}
