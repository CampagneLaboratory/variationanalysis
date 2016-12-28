package org.campagnelab.dl.framework.mappers;

import org.nd4j.linalg.api.ndarray.INDArray;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.util.HashMap;
import java.util.function.Function;

/**
 * Maps a list of labels to a fixed set of integers, which are used to encode each position in one-hot encoding,
 * at a specific baseIndex into the list of labels. Labels are assumed to range from 0 to numLabels - 1.
 * Created by joshuacohen on 11/21/16.
 */
public class OneHotBaseLabelMapper<RecordType> implements LabelMapper<RecordType> {
    static private Logger LOG = LoggerFactory.getLogger(OneHotBaseFeatureMapper.class);

    private int baseIndex;
    private int numLabels;
    private Function<RecordType, int[]> recordToLabel;
    private HashMap<RecordType, int[]> recordLabelMap;

    private static final int[] indices = new int[]{0, 0};

    /**
     * Creates a OneHotBaseLabelMapper with a specified baseIndex and conversion function
     * @param baseIndex the base index into the set of labels that this mapper is mapping
     * @param numLabels the number of labels that the record can correspond to
     * @param recordToLabel function that takes in a record and outputs an int array, corresponding to
     *                      the labels (represented by numbers ranging from 0 to numLabels - 1) at each
     *                      position in the record
     */
    public OneHotBaseLabelMapper(int baseIndex, int numLabels, Function<RecordType, int[]> recordToLabel) {
        this.baseIndex = baseIndex;
        this.numLabels = numLabels;
        this.recordToLabel = recordToLabel;
        recordLabelMap = new HashMap<>();
    }

    @Override
    public void prepareToNormalize(RecordType record, int indexOfRecord) {
        recordLabelMap.put(record, recordToLabel.apply(record));
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
        for (int labelIndex = 0; labelIndex < numberOfLabels(); labelIndex++) {
            indices[1] = labelIndex;
            labels.putScalar(indices, produceLabel(record, labelIndex));
        }
    }

    @Override
    public float produceLabel(RecordType record, int labelIndex) {
        int[] label = recordLabelMap.computeIfAbsent(record, r -> {
            int[] l = recordToLabel.apply(r);
            recordLabelMap.put(r, l);
            return l;
        });
        if (baseIndex < 0 || baseIndex >= label.length) {
            LOG.warn("incompatible base index: {} for label: {} of length {}",
                    baseIndex, label, label.length);
            return numLabels - 1;
        }
        return labelIndex == label[baseIndex] ? 1F : 0F;
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
