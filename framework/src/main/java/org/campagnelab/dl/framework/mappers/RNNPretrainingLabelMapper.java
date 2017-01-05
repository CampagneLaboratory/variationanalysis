package org.campagnelab.dl.framework.mappers;

import org.nd4j.linalg.api.ndarray.INDArray;

/**
 * Deprecated implementation of a label mapper for pretraining- used to compare to implementation
 * with concatenating mappers
 *
 * Maps a feature sequence ACT to the corresponding sequence "...ACT!", where ! is the EOS token
 * and . is a blank space token, that is masked 0
 *
 * Created by joshuacohen on 12/8/16.
 */
@Deprecated
public class RNNPretrainingLabelMapper<RecordType> implements LabelMapper<RecordType> {
    private RNNFeatureMapper<RecordType> labelMapper;
    private MappedDimensions dim;
    private int labelsPerTimeStep;
    private int maxSequenceLen;

    private int[] mapperIndices = new int[]{0, 0, 0};
    private int[] maskerIndices = new int[]{0, 0};

    public RNNPretrainingLabelMapper(RNNFeatureMapper<RecordType> labelMapper) {
        this.labelMapper = labelMapper;
        labelsPerTimeStep = labelMapper.featuresPerTimeStep();
        maxSequenceLen = labelMapper.maxSequenceLength();
        dim = new MappedDimensions(labelsPerTimeStep,
                maxSequenceLen * 2 + 1);
    }

    @Override
    public int numberOfLabels() {
        return labelMapper.numberOfFeatures() * 2 + 1;
    }

    @Override
    public MappedDimensions dimensions() {
        return dim;
    }

    @Override
    public void prepareToNormalize(RecordType record, int indexOfRecord) {
        labelMapper.prepareToNormalize(record, indexOfRecord);
    }

    @Override
    public void mapLabels(RecordType record, INDArray inputs, int indexOfRecord) {
        labelMapper.prepareToNormalize(record, indexOfRecord);
        mapperIndices[0] = indexOfRecord;
        for (int i = 0; i < maxSequenceLen; i++) {
            mapperIndices[2] = i;
            int indexInSequence = i % labelMapper.sequenceLength;
            for (int j = 0; j < labelsPerTimeStep; j++) {
                mapperIndices[1] = j;
                int labelIndex = indexInSequence * labelsPerTimeStep + j;
                inputs.putScalar(mapperIndices, produceLabel(record, labelIndex));
            }
        }
    }

    @Override
    public boolean hasMask() {
        return true;
    }

    @Override
    public void maskLabels(RecordType record, INDArray mask, int indexOfRecord) {
        labelMapper.prepareToNormalize(record, indexOfRecord);
        maskerIndices[0] = indexOfRecord;
        int sequenceLenWithPadding = sequenceLenWithPadding(labelMapper.sequenceLength);
        for (int i = 0; i < maxSequenceLen; i++) {
            maskerIndices[1] = i;
            // TODO : Check offsets
            mask.putScalar(maskerIndices, i < labelMapper.sequenceLength || i >= sequenceLenWithPadding ? 0F : 1F);
        }
    }

    @Override
    public boolean isMasked(RecordType record, int labelIndex) {
        int sequenceLenWithPadding = sequenceLenWithPadding(labelMapper.sequenceLength);
        int timeStepIdx = labelIndex / labelsPerTimeStep;
        // TODO : Check offsets
        return timeStepIdx >= labelMapper.sequenceLength && timeStepIdx < sequenceLenWithPadding;
    }

    @Override
    public float produceLabel(RecordType record, int labelIndex) {
        // TODO : Check offsets
        int sequenceLenWithPadding = sequenceLenWithPadding(labelMapper.sequenceLength);
        int timeStepIdx = labelIndex / labelsPerTimeStep;
        return (timeStepIdx >= labelMapper.sequenceLength && timeStepIdx < sequenceLenWithPadding)
                ? labelMapper.produceFeature(record, labelIndex)
                : 0F;
    }

    private int sequenceLenWithPadding(int sequenceLen) {
        return sequenceLen * 2;
    }
}
