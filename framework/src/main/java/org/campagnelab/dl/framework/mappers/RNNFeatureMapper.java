package org.campagnelab.dl.framework.mappers;

import org.nd4j.linalg.api.ndarray.INDArray;

import java.util.function.Function;

/**
 * Concatenating mapper that takes in either a set of delegate mappers that map to 1 dimensional
 * features, or constructs a set of OneHotBaseFeatureMappers to use, and maps to two dimensional
 * features that can be used in recurrent neural networks. Delegates are assumed to map for
 * one base into the sequence. The number of features for each delegate corresponds to the
 * number of features per time step in the two-dimensional features produced by the mapper,
 * and the number of delegates- or the maxSequenceLength, which creates a set of
 * OneHotBaseFeatureMappers with cardinality = maxSequenceLength- is the number of time
 * steps per feature.
 *
 * For sequences with lengths that are less than the number of delegates (i.e., the maxSequenceLength),
 * extra base indexes will not be masked- i.e., for # delegates = 6, and a sequence "AC", any features corresponding
 * to the 2nd through 5th indices will be 0, with a mask of 0
 *
 * Example-
 * Sequence = ACT
 * OneHotBaseMapper at baseIndex = 0- maps to [1, 0, 0, 0]
 * OneHotBaseMapper at baseIndex = 1- maps to [0, 1, 0, 0]
 * OneHotBaseMapper at baseIndex = 2- maps to [0, 0, 1, 0]
 * RNNFeatureMapper-
 * [[1, 0, 0],
 *  [0, 1, 0],
 *  [0, 0, 1],
 *  [0, 0, 0]]
 * Created by joshuacohen on 11/18/16.
 */
public class RNNFeatureMapper<RecordType> implements FeatureMapper<RecordType> {
    private int featuresPerTimeStep;
    private FeatureMapper<RecordType>[] delegates;
    private Function<RecordType, Integer> recordToSequenceLength;

    private int[] indicesMapper = new int[]{0, 0, 0};
    private int[] indicesMasker = new int[]{0, 0};

    int sequenceLength;

    private MappedDimensions dim;

    /**
     * Constructor used to create an RNNFeatureMapper with a set of OneHotBaseFeatureMappers with
     * cardinality = maxSequenceLength as the delegate 1-dimensional feature mappers
     * @param maxSequenceLength Maximum sequence length of any sequence in the set
     * @param recordToString Function that converts a record to a string representation- used
     *                       in OneHotBaseFeatureMapper
     * @param recordToSequenceLength Function that converts a record of a sequence to its length
     */
    public RNNFeatureMapper(int maxSequenceLength, Function<RecordType, String> recordToString,
                            Function<RecordType, Integer> recordToSequenceLength) {
        this(recordToSequenceLength,
                createOneHotBaseFeatureMappers(maxSequenceLength, recordToString));
    }

    /**
     * Constructor used to create an RNNFeatureMapper with an arbitrary set of delegates, all
     * with the same number of features
     * @param recordToSequenceLength Function that converts a record of a sequence to its length
     * @param delegates Delegate feature mappers
     */
    @SafeVarargs
    public RNNFeatureMapper(Function<RecordType, Integer> recordToSequenceLength,
                            FeatureMapper<RecordType>... delegates) {
        this.recordToSequenceLength = recordToSequenceLength;
        MappedDimensions dimensions = delegates[0].dimensions();
        for (FeatureMapper<RecordType> featureMapper : delegates) {
            MappedDimensions featureMapperDimensions = featureMapper.dimensions();
            if (!dimensions.equals(featureMapperDimensions)) {
                throw new RuntimeException("All delegate feature mappers should have same dimensions");
            }
            if (featureMapperDimensions.numDimensions() != 1) {
                throw new RuntimeException("All delegate feature mappers should be one dimensional");
            }
        }
        this.featuresPerTimeStep = dimensions.numElements();
        this.delegates = delegates;
        dim = new MappedDimensions(featuresPerTimeStep, delegates.length);
    }

    @Override
    public int numberOfFeatures() {
        return delegates.length * featuresPerTimeStep;
    }

    @Override
    public MappedDimensions dimensions() {
        assert dim.numElements() == numberOfFeatures() : "Number of elements must match number of features.";
        return dim;
    }

    @Override
    public void prepareToNormalize(RecordType record, int indexOfRecord) {
        sequenceLength = recordToSequenceLength.apply(record);
        for (FeatureMapper<RecordType> delegate : delegates) {
            delegate.prepareToNormalize(record, indexOfRecord);
        }
    }

    @Override
    public void mapFeatures(RecordType record, INDArray inputs, int indexOfRecord) {
        indicesMapper[0] = indexOfRecord;
        for (int i = 0; i < delegates.length; i++) {
            indicesMapper[2] = i;
            final FeatureMapper<RecordType> delegate = delegates[i];
            for (int j = 0; j < delegate.numberOfFeatures(); j++) {
                indicesMapper[1] = j;
                if (i < sequenceLength) {
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
            mask.putScalar(indicesMasker, i < sequenceLength ? 1F : 0F);
        }
    }

    @Override
    public float produceFeature(RecordType record, int featureIndex) {
        int delegateIdx = featureIndex / featuresPerTimeStep;
        if (delegateIdx >= sequenceLength) {
            return 0F;
        }
        int featureInDelegateIdx = featureIndex % featuresPerTimeStep;
        return delegates[delegateIdx].produceFeature(record, featureInDelegateIdx);
    }

    @Override
    public boolean isMasked(RecordType record, int featureIndex) {
        int delegateIdx = featureIndex / featuresPerTimeStep;
        if (delegateIdx >= sequenceLength) {
            return false;
        }
        int featureInDelegateIdx = featureIndex % featuresPerTimeStep;
        return !delegates[delegateIdx].hasMask()
                || delegates[delegateIdx].isMasked(record, featureInDelegateIdx);
    }

    int maxSequenceLength() {
        return delegates.length;
    }

    int featuresPerTimeStep() {
        return featuresPerTimeStep;
    }

    private static <RecordType> FeatureMapper<RecordType>[] createOneHotBaseFeatureMappers(int maxSequenceLength,
                                                           Function<RecordType, String> recordToString) {
        FeatureMapper<RecordType>[] delegates = new OneHotBaseFeatureMapper[maxSequenceLength];
        for (int i = 0; i < maxSequenceLength; i++) {
            delegates[i] = new OneHotBaseFeatureMapper<>(i, recordToString);
        }
        return delegates;
    }
}
