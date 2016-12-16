package org.campagnelab.dl.framework.mappers;

import org.campagnelab.dl.framework.mappers.FeatureMapper;
import org.campagnelab.dl.framework.mappers.MappedDimensions;
import org.nd4j.linalg.api.ndarray.INDArray;

import java.lang.reflect.Array;
import java.util.Arrays;

/**
 * Concatenates 2D features from a set of feature mappers horizontally. Delegate feature mappers
 * must have the same dimensionality of 2 and have the same number of features- i.e., the same
 * number of features per time step, corresponding to the first dimension in the 2D features.
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
 * The TwoDimensionalConcatFeatureMapper will map the following
 * [[1, 0, 0],
 *  [0, 1, 0],
 *  [0, 0, 1]]
 *
 * Created by joshuacohen on 12/8/16.
 */
public class TwoDimensionalConcatFeatureMapper<RecordType> implements FeatureMapper<RecordType> {
    private FeatureMapper<RecordType>[] delegates;
    private int featuresPerTimeStep;
    private int numFeatures;
    private int[] timeStepOffsets;
    private int[] numFeaturesOffsets;
    private MappedDimensions dim;
    private int[] mapperIndices = new int[]{0, 0, 0};
    private int[] maskerIndices = new int[]{0, 0};
    private int zeroPaddingWidth;
    private int totalTimeSteps;

    /**
     * Creates a concatenating feature mapper from a set of delegate mappers with padding zeros ended
     * to the end of each feature
     * @param zeroPaddingWidth Amount of zeros to add to end of feature
     * @param delegates        Feature mappers that are concatenated
     */
    @SafeVarargs
    public TwoDimensionalConcatFeatureMapper(int zeroPaddingWidth, FeatureMapper<RecordType>... delegates) {
        this.zeroPaddingWidth = zeroPaddingWidth;
        MappedDimensions dimensions = delegates[0].dimensions();
        timeStepOffsets = new int[delegates.length + 1];
        numFeaturesOffsets = new int[delegates.length + 1];
        timeStepOffsets[0] = 0;
        numFeaturesOffsets[0] = 0;
        int totalTimeSteps = 0;
        numFeatures = 0;
        int i = 1;
        for (FeatureMapper<RecordType> delegate : delegates) {
            MappedDimensions delegateDimensions = delegate.dimensions();
            if (delegateDimensions.numDimensions() != 2) {
                throw new RuntimeException("Delegate mappers must be two dimensional mappers");
            }
            if (!dimensions.equalsDimension(delegateDimensions, 1)) {
                throw new RuntimeException("Delegate mappers must have same number of features");
            }
            int delegateTimeSteps = delegate.dimensions().numElements(2);
            numFeatures += delegate.numberOfFeatures();
            totalTimeSteps += delegateTimeSteps;
            timeStepOffsets[i] = delegateTimeSteps;
            numFeaturesOffsets[i] = numFeatures;
            i++;
        }
        this.totalTimeSteps = totalTimeSteps;
        featuresPerTimeStep = dimensions.numElements(1);
        this.delegates = delegates;
        dim = new MappedDimensions(featuresPerTimeStep + zeroPaddingWidth, totalTimeSteps);
    }

    /**
     * Creates a concatenating feature mapper from a set of delegate mappers with no padding
     * @param delegates Feature mappers that are concatenated
     */
    @SafeVarargs
    public TwoDimensionalConcatFeatureMapper(FeatureMapper<RecordType>... delegates) {
        this(0, delegates);
    }

    @Override
    public int numberOfFeatures() {
        return numFeatures + (totalTimeSteps * zeroPaddingWidth);
    }

    @Override
    public MappedDimensions dimensions() {
        return dim;
    }

    @Override
    public void prepareToNormalize(RecordType record, int indexOfRecord) {

    }

    @Override
    public void mapFeatures(RecordType record, INDArray inputs, int indexOfRecord) {
        mapperIndices[0] = indexOfRecord;
        int offset = 0;
        for (FeatureMapper<RecordType> delegate : delegates) {
            int numTimeSteps = delegate.dimensions().numElements(2);
            for (int i = 0; i < numTimeSteps; i++) {
                mapperIndices[2] = offset;
                for (int j = 0; j < featuresPerTimeStep + zeroPaddingWidth; j++) {
                    mapperIndices[1] = j;
                    int featureIndex = i * featuresPerTimeStep + j;
                    float feature = j < featuresPerTimeStep ? delegate.produceFeature(record, featureIndex) : 0F;
                    inputs.putScalar(mapperIndices, feature);
                }
                offset++;
            }
        }
    }

    @Override
    public boolean hasMask() {
        return true;
    }

    @Override
    public void maskFeatures(RecordType record, INDArray mask, int indexOfRecord) {
        maskerIndices[0] = indexOfRecord;
        int offset = 0;
        for (FeatureMapper<RecordType> delegate : delegates) {
            int numTimeSteps = delegate.dimensions().numElements(2);
            for (int i = 0; i < numTimeSteps; i++) {
                maskerIndices[1] = offset;
                mask.putScalar(maskerIndices, delegate.isMasked(record, i * featuresPerTimeStep) ? 1F : 0F);
                offset++;
            }
        }
    }

    @Override
    public boolean isMasked(RecordType record, int featureIndex) {
        int timeStepIndex = featureIndex / (featuresPerTimeStep + zeroPaddingWidth);
        int featureInTimeStepIndex = featureIndex % (featuresPerTimeStep + zeroPaddingWidth);
        int newFeatureIndex = timeStepIndex * featuresPerTimeStep + featureInTimeStepIndex;
        int indexOfDelegate = Arrays.binarySearch(timeStepOffsets, timeStepIndex);
        if (indexOfDelegate < 0) {
            indexOfDelegate = -(indexOfDelegate + 1) - 1;
        }
        return delegates[indexOfDelegate].isMasked(record,
                newFeatureIndex - numFeaturesOffsets[indexOfDelegate]);
    }

    @Override
    public float produceFeature(RecordType record, int featureIndex) {
        int timeStepIndex = featureIndex / (featuresPerTimeStep + zeroPaddingWidth);
        int featureInTimeStepIndex = featureIndex % (featuresPerTimeStep + zeroPaddingWidth);
        int newFeatureIndex = timeStepIndex * featuresPerTimeStep + featureInTimeStepIndex;
        int indexOfDelegate = Arrays.binarySearch(timeStepOffsets, newFeatureIndex);
        if (indexOfDelegate < 0) {
            indexOfDelegate = -(indexOfDelegate + 1) - 1;
        }
        boolean validFeature = featureInTimeStepIndex < featuresPerTimeStep;
        try {
            Object a = delegates[indexOfDelegate];
            int b = timeStepOffsets[indexOfDelegate];
        } catch (ArrayIndexOutOfBoundsException e) {
            int c = 1;
        }
        return validFeature ? delegates[indexOfDelegate].produceFeature(record,
                newFeatureIndex - numFeaturesOffsets[indexOfDelegate]) : 0F;
    }
}
