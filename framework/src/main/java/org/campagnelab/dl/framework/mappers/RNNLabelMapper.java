package org.campagnelab.dl.framework.mappers;

import org.nd4j.linalg.api.ndarray.INDArray;

import java.util.function.Function;

/**
 * Concatenating mapper that takes in either a set of delegate mappers that map to 1 dimensional
 * labels, or constructs a set of OneHotBaseLabelMappers to use, and maps to two dimensional
 * labels that can be used in recurrent neural networks. Delegates are assumed to map for
 * one base into the sequence. The number of labels for each delegate corresponds to the
 * number of labels per time step in the two-dimensional features produced by the mapper,
 * and the number of delegates- or the maxSequenceLength, which creates a set of
 * OneHotBaseLabelMappers with cardinality = maxSequenceLength- is the number of time
 * steps per label.
 *
 * For sequences with lengths that are less than the number of delegates (i.e., the maxSequenceLength),
 * extra base indexes will not be masked- i.e., for # delegates = 6, and a label sequence [0,1], any features
 * corresponding to the 2nd through 5th indices will be 0, with a mask of 0
 *
 * Example-
 * Sequence = ACT with a label [0, 1, 2]
 * OneHotBaseMapper at baseIndex = 0- maps to [1, 0, 0, 0]
 * OneHotBaseMapper at baseIndex = 1- maps to [0, 1, 0, 0]
 * OneHotBaseMapper at baseIndex = 2- maps to [0, 0, 1, 0]
 * RNNFeatureMapper-
 * [[1, 0, 0],
 *  [0, 1, 0],
 *  [0, 0, 1],
 *  [0, 0, 0]]
 * Created by joshuacohen on 11/21/16.
 */
public class RNNLabelMapper<RecordType> implements LabelMapper<RecordType> {
    private int labelsPerTimeStep;
    private LabelMapper<RecordType>[] delegates;
    Function<RecordType, Integer> recordToSequenceLength;

    private int[] indicesMapper = new int[]{0, 0, 0};
    private int[] indicesMasker = new int[]{0, 0};

    private MappedDimensions dim;

    private int sequenceLength;

    /**
     * Constructor used to create an RNNLabelMapper with a set of OneHotBaseLabelMappers with
     * cardinality = maxSequenceLength as the delegate 1-dimensional feature mappers
     * @param maxSequenceLength Maximum sequence length of any sequence in the set
     * @param labelsPerTimeStep Number of labels for each OneHotBaseLabelMapper
     * @param recordToLabel Function that converts a record to a string representation- used
     *                       in OneHotBaseLabelMapper
     * @param recordToSequenceLength Function that converts a record of a sequence to its length
     */
    public RNNLabelMapper(int maxSequenceLength, int labelsPerTimeStep, Function<RecordType, int[]> recordToLabel,
                            Function<RecordType, Integer> recordToSequenceLength) {
        this(recordToSequenceLength,
                createOneHotBaseLabelMappers(maxSequenceLength, labelsPerTimeStep, recordToLabel));
    }

    /**
     * Constructor used to create an RNNFeatureMapper with an arbitrary set of delegates, all
     * with the same number of features
     * @param recordToSequenceLength Function that converts a record of a sequence to its length
     * @param delegates Delegate label mappers
     */
    @SafeVarargs
    public RNNLabelMapper(Function<RecordType, Integer> recordToSequenceLength,
                          LabelMapper<RecordType>... delegates) {
        this.recordToSequenceLength = recordToSequenceLength;
        MappedDimensions dimensions = delegates[0].dimensions();
        for (LabelMapper<RecordType> labelMapper : delegates) {
            MappedDimensions labelMapperDimensions = labelMapper.dimensions();
            if (!dimensions.equals(labelMapperDimensions)) {
                throw new RuntimeException("All delegate label mappers should have same dimensions");
            }
            if (labelMapperDimensions.numDimensions() != 1) {
                throw new RuntimeException("All delegate label mappers should be one dimensional");
            }
        }
        this.labelsPerTimeStep = dimensions.numElements();
        this.delegates = delegates;
        dim = new MappedDimensions(labelsPerTimeStep, delegates.length);
    }

    @Override
    public int numberOfLabels() {
        return delegates.length * labelsPerTimeStep;
    }

    @Override
    public MappedDimensions dimensions() {
        assert dim.numElements() == numberOfLabels() : "Number of elements must match number of labels.";
        return dim;
    }

    @Override
    public void mapLabels(RecordType record, INDArray labels, int indexOfRecord) {
        indicesMapper[0] = indexOfRecord;
        for (int i = 0; i < delegates.length; i++) {
            indicesMapper[2] = i;
            final LabelMapper<RecordType> delegate = delegates[i];
            for (int j = 0; j < delegates[i].numberOfLabels(); j++) {
                indicesMapper[1] = j;
                if (i < sequenceLength) {
                    labels.putScalar(indicesMapper, delegate.produceLabel(record, j));
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
            mask.putScalar(indicesMasker, i < sequenceLength ? 1F : 0F);
        }
    }

    @Override
    public float produceLabel(RecordType record, int labelIndex) {
        int delegateIdx = labelIndex / labelsPerTimeStep;
        if (delegateIdx >= sequenceLength) {
            return 0F;
        }
        int labelInDelegateIdx = labelIndex % labelsPerTimeStep;
        return delegates[delegateIdx].produceLabel(record, labelInDelegateIdx);
    }

    @Override
    public boolean isMasked(RecordType record, int labelIndex) {
        int delegateIdx = labelIndex / labelsPerTimeStep;
        if (delegateIdx >= sequenceLength) {
            return false;
        }
        int labelInDelegateIdx = labelIndex % labelsPerTimeStep;
        return !delegates[delegateIdx].hasMask()
                || delegates[delegateIdx].isMasked(record, labelInDelegateIdx);
    }

    @Override
    public void prepareToNormalize(RecordType record, int indexOfRecord) {
        sequenceLength = recordToSequenceLength.apply(record);
        for (LabelMapper<RecordType> delegate : delegates) {
            delegate.prepareToNormalize(record, indexOfRecord);
        }
    }

    public int maxSequenceLength() {
        return delegates.length;
    }

    public int labelsPerTimeStep() {
        return labelsPerTimeStep;
    }

    private static <RecordType> LabelMapper<RecordType>[] createOneHotBaseLabelMappers(int maxSequenceLength,
                                                                                       int numLabels,
                                                                                       Function<RecordType, int[]> recordToLabel) {
        LabelMapper<RecordType>[] delegates = new OneHotBaseLabelMapper[maxSequenceLength];
        for (int i = 0; i < maxSequenceLength; i++) {
            delegates[i] = new OneHotBaseLabelMapper<>(i, numLabels, recordToLabel);
        }
        return delegates;
    }
}
