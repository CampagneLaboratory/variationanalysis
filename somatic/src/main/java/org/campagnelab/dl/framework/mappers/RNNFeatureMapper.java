package org.campagnelab.dl.framework.mappers;

import org.nd4j.linalg.api.ndarray.INDArray;

import java.util.function.Function;

/**
 * Created by joshuacohen on 11/18/16.
 */
public class RNNFeatureMapper<RecordType> implements FeatureMapper<RecordType> {
    private int featuresPerTimeStep;
    private OneHotBaseFeatureMapper<RecordType>[] delegates;
    private Function<RecordType, Integer> recordToSequenceLength;

    private int[] indicesMapper = new int[]{0, 0, 0};
    private int[] indicesMasker = new int[]{0, 0};

    public RNNFeatureMapper(int featuresPerTimeStep, int maxSequenceLength, Function<RecordType, String> recordToString,
                            Function<RecordType, Integer> recordToSequenceLength) {
        this.recordToSequenceLength = recordToSequenceLength;
        delegates = new OneHotBaseFeatureMapper[maxSequenceLength];
        for (int i = 0; i < maxSequenceLength; i++) {
            delegates[i] = new OneHotBaseFeatureMapper<>(i, recordToString);
        }
        this.featuresPerTimeStep = delegates[0].numberOfFeatures();
        for (OneHotBaseFeatureMapper<RecordType> oneHotBaseFeatureMapper : delegates) {
            if (oneHotBaseFeatureMapper.numberOfFeatures() != featuresPerTimeStep) {
                throw new RuntimeException("All delegate one hot base mappers should have same number of features");
            }
        }
    }

    @Override
    public int numberOfFeatures() {
        return delegates.length * featuresPerTimeStep;
    }

    @Override
    public MappedDimensions dimensions() {

        // TODO Josh, I think the next line should use delegates.length. not numberOfFeatures().
        MappedDimensions dim = new MappedDimensions(featuresPerTimeStep, numberOfFeatures());
        assert dim.numElements() == numberOfFeatures() : "Number of elements must match number of features.";
        // TODO: build dim in the constructor and return the same object every time. It is a constant.
        return dim;
    }

    @Override
    public void prepareToNormalize(RecordType record, int indexOfRecord) {

    }

    @Override
    public void mapFeatures(RecordType record, INDArray inputs, int indexOfRecord) {
        indicesMapper[0] = indexOfRecord;
        for (int i = 0; i < delegates.length; i++) {
            indicesMapper[2] = i;
            final OneHotBaseFeatureMapper<RecordType> delegate = delegates[i];
            for (int j = 0; j < delegate.numberOfFeatures(); j++) {
                indicesMapper[1] = j;
                if (i < recordToSequenceLength.apply(record)) {
                    inputs.putScalar(indicesMapper, delegate.produceFeature(record, j));
                } else {
                    inputs.putScalar(indicesMapper, 0F);
                }
            }
        }
    }

    @Override
    public boolean hasMask() {
        return true;
    }

    @Override
    public void maskFeatures(RecordType record, INDArray mask, int indexOfRecord) {
        indicesMasker[0] = indexOfRecord;
        for (int i = 0; i < delegates.length; i++) {
            indicesMasker[1] = i;
            mask.putScalar(indicesMasker, i < recordToSequenceLength.apply(record) ? 1F : 0F);
        }
    }

    @Override
    public float produceFeature(RecordType record, int featureIndex) {
        int delegateIdx = featureIndex / featuresPerTimeStep;
        int featureInDelegateIdx = featureIndex % featuresPerTimeStep;
        return delegates[delegateIdx].produceFeature(record, featureInDelegateIdx);
    }

    @Override
    public boolean isMasked(RecordType record, int featureIndex) {
        int delegateIdx = featureIndex / featuresPerTimeStep;
        int featureInDelegateIdx = featureIndex % featuresPerTimeStep;
        return delegates[delegateIdx].isMasked(record, featureInDelegateIdx);
    }
}
