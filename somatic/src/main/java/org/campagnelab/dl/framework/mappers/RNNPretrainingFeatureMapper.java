package org.campagnelab.dl.framework.mappers;

import org.nd4j.linalg.api.ndarray.INDArray;

/**
 * Deprecated implementation of a feature mapper for pretraining- used to compare to implementation
 * with concatenating mappers
 *
 * Maps a feature sequence ACT to the corresponding sequence "ACT!ACT", where ! is the EOS token
 *
 * Created by joshuacohen on 12/7/16.
 */
@Deprecated
public class RNNPretrainingFeatureMapper<RecordType> implements FeatureMapper<RecordType> {
    private RNNFeatureMapper<RecordType> featureMapper;
    private MappedDimensions dim;
    private int featuresPerTimeStep;
    private int maxSequenceLen;

    private int[] mapperIndices = new int[]{0, 0, 0};
    private int[] maskerIndices = new int[]{0, 0};

    public RNNPretrainingFeatureMapper(RNNFeatureMapper<RecordType> featureMapper) {
        this.featureMapper = featureMapper;
        featuresPerTimeStep = featureMapper.featuresPerTimeStep();
        maxSequenceLen = featureMapper.maxSequenceLength();
        dim = new MappedDimensions(featuresPerTimeStep,
                maxSequenceLen * 2 + 1);
    }

    @Override
    public int numberOfFeatures() {
        return featureMapper.numberOfFeatures() * 2 + 1;
    }

    @Override
    public MappedDimensions dimensions() {
        return dim;
    }

    @Override
    public void prepareToNormalize(RecordType record, int indexOfRecord) {
        featureMapper.prepareToNormalize(record, indexOfRecord);
    }

    @Override
    public void mapFeatures(RecordType record, INDArray inputs, int indexOfRecord) {
        featureMapper.prepareToNormalize(record, indexOfRecord);
        mapperIndices[0] = indexOfRecord;
        int sequenceLen = featureMapper.sequenceLengths[indexOfRecord];
        for (int i = 0; i < maxSequenceLen; i++) {
            mapperIndices[2] = i;
            int indexInSequence = i % sequenceLen;
            for (int j = 0; j < featuresPerTimeStep; j++) {
                mapperIndices[1] = j;
                int featureIndex = indexInSequence * featuresPerTimeStep + j;
                inputs.putScalar(mapperIndices, produceFeature(record, featureIndex));
            }
        }
    }

    @Override
    public boolean hasMask() {
        return true;
    }

    @Override
    public void maskFeatures(RecordType record, INDArray mask, int indexOfRecord) {
        featureMapper.prepareToNormalize(record, indexOfRecord);
        maskerIndices[0] = indexOfRecord;
        int sequenceLenWithPadding = sequenceLenWithPadding(featureMapper.sequenceLengths[indexOfRecord]);
        for (int i = 0; i < maxSequenceLen; i++) {
            maskerIndices[1] = i;
            // TODO: Handling of EOS?
            mask.putScalar(maskerIndices, i < sequenceLenWithPadding ? 1F : 0F);
        }
    }

    @Override
    public boolean isMasked(RecordType record, int featureIndex) {
        Integer sequenceLengthFromMap = featureMapper.sequenceLengthMap.get(record);
        int sequenceLength = sequenceLengthFromMap != null
                ? sequenceLengthFromMap : featureMapper.recordToSequenceLength.apply(record);
        int sequenceLenWithPadding = sequenceLenWithPadding(sequenceLength);
        int timeStepIdx = featureIndex / featuresPerTimeStep;
        return timeStepIdx < sequenceLenWithPadding;
    }

    @Override
    public float produceFeature(RecordType record, int featureIndex) {
        Integer sequenceLengthFromMap = featureMapper.sequenceLengthMap.get(record);
        int sequenceLength = sequenceLengthFromMap != null
                ? sequenceLengthFromMap : featureMapper.recordToSequenceLength.apply(record);
        int sequenceLenWithPadding = sequenceLenWithPadding(sequenceLength);
        int timeStepIdx = featureIndex / featuresPerTimeStep;
        // TODO : Handling of EOS?
        return (timeStepIdx < sequenceLenWithPadding) ? featureMapper.produceFeature(record, featureIndex)
                : 0F;
    }

    private int sequenceLenWithPadding(int sequenceLen) {
        return sequenceLen * 2 + 1;
    }
}
