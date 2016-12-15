package org.campagnelab.dl.framework.mappers;

import org.campagnelab.dl.framework.mappers.LabelMapper;
import org.campagnelab.dl.framework.mappers.MappedDimensions;
import org.nd4j.linalg.api.ndarray.INDArray;

import java.util.Arrays;

/**
 * Concatenates 2D labels from a set of label mappers horizontally. Delegate label mappers
 * must have the same dimensionality of 2 and have the same number of labels- i.e., the same
 * number of labels per time step, corresponding to the first dimension in the 2D labels.
 *
 * For example, if there are three mappers, that map as follows
 *
 * Mapper 1-
 * [[1],
 *  [0],
 *  [0]]
 *
 * Mapper 2-
 * [[0],
 *  [1],
 *  [0]]
 *
 * Mapper 3-
 * [[0],
 *  [0],
 *  [1]]
 *
 * The TwoDimensionalConcatLabelMapper will map the following
 * [[1, 0, 0],
 *  [0, 1, 0],
 *  [0, 0, 1]]
 *
 * Created by joshuacohen on 12/8/16.
 */
public class TwoDimensionalConcatLabelMapper<RecordType> implements LabelMapper<RecordType> {
    private LabelMapper<RecordType>[] delegates;
    private int labelsPerTimeStep;
    private int numLabels;
    private int[] numLabelsOffsets;
    private MappedDimensions dim;
    private int[] mapperIndices = new int[]{0, 0, 0};
    private int[] maskerIndices = new int[]{0, 0};
    private int zeroPaddingWidth;


    /**
     * Creates a concatenating label mapper from a set of delegate mappers with padding zeros added
     * to the end of each label
     * @param zeroPaddingWidth Amount of zeros to add to end of label
     * @param delegates        Label mappers that are concatenated
     */
    @SafeVarargs
    public TwoDimensionalConcatLabelMapper(int zeroPaddingWidth, LabelMapper<RecordType>... delegates) {
        this.zeroPaddingWidth = zeroPaddingWidth;
        MappedDimensions dimensions = delegates[0].dimensions();
        numLabelsOffsets = new int[delegates.length + 1];
        numLabelsOffsets[0] = 0;
        int totalTimeSteps = 0;
        numLabels = 0;
        int i = 1;
        for (LabelMapper<RecordType> delegate : delegates) {
            MappedDimensions delegateDimensions = delegate.dimensions();
            if (delegateDimensions.numDimensions() != 2) {
                throw new RuntimeException("Delegate mappers must be two dimensional mappers");
            }
            if (!dimensions.equalsDimension(delegateDimensions, 1)) {
                throw new RuntimeException("Delegate mappers must have same number of labels");
            }
            int delegateTimeSteps = delegate.dimensions().numElements(2);
            numLabels += delegate.numberOfLabels() + delegateTimeSteps * zeroPaddingWidth;
            totalTimeSteps += delegateTimeSteps;
            numLabelsOffsets[i] = numLabels;
            i++;
        }
        labelsPerTimeStep = dimensions.numElements(1);
        this.delegates = delegates;
        dim = new MappedDimensions(labelsPerTimeStep + zeroPaddingWidth, totalTimeSteps);
    }

    /**
     * Creates a concatenating label mapper from a set of delegate mappers with no padding
     * @param delegates Label mappers that are concatenated
     */
    @SafeVarargs
    public TwoDimensionalConcatLabelMapper(LabelMapper<RecordType>... delegates) {
        this(0, delegates);
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
    public void mapLabels(RecordType record, INDArray labels, int indexOfRecord) {
        mapperIndices[0] = indexOfRecord;
        int offset = 0;
        for (LabelMapper<RecordType> delegate : delegates) {
            int numTimeSteps = delegate.dimensions().numElements(2);
            for (int i = 0; i < numTimeSteps; i++) {
                mapperIndices[2] = offset;
                for (int j = 0; j < labelsPerTimeStep + zeroPaddingWidth; j++) {
                    mapperIndices[1] = j;
                    int featureIdx = i * labelsPerTimeStep + j;
                    float label = j <= labelsPerTimeStep ? delegate.produceLabel(record, featureIdx) : 0F;
                    labels.putScalar(mapperIndices, label);
                }
                offset++;
            }
        }
    }

    @Override
    public float produceLabel(RecordType record, int labelIndex) {
        int indexOfDelegate = Arrays.binarySearch(numLabelsOffsets, labelIndex);
        if (indexOfDelegate < 0) {
            indexOfDelegate = -(indexOfDelegate + 1) - 1;
        }
        float label = delegates[indexOfDelegate].produceLabel(record,
                labelIndex - numLabelsOffsets[indexOfDelegate]);
        boolean labelValid = labelIndex % (labelsPerTimeStep + zeroPaddingWidth) <= labelsPerTimeStep;
        return labelValid ? label : 0F;
    }

    @Override
    public boolean hasMask() {
        return true;
    }

    @Override
    public void maskLabels(RecordType record, INDArray mask, int indexOfRecord) {
        maskerIndices[0] = indexOfRecord;
        int offset = 0;
        for (LabelMapper<RecordType> delegate : delegates) {
            int numTimeSteps = delegate.dimensions().numElements(2);
            for (int i = 0; i < numTimeSteps; i++) {
                maskerIndices[1] = offset;
                mask.putScalar(maskerIndices, delegate.isMasked(record, i * labelsPerTimeStep) ? 1F : 0F);
                offset++;
            }
        }
    }

    @Override
    public boolean isMasked(RecordType record, int featureIndex) {
        int indexOfDelegate = Arrays.binarySearch(numLabelsOffsets, featureIndex);
        if (indexOfDelegate < 0) {
            indexOfDelegate = -(indexOfDelegate + 1) - 1;
        }
        return delegates[indexOfDelegate].isMasked(record,
                featureIndex - numLabelsOffsets[indexOfDelegate]);
    }

    @Override
    public void prepareToNormalize(RecordType record, int indexOfRecord) {

    }
}
