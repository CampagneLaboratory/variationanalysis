package org.campagnelab.dl.framework.mappers.processing;

import org.campagnelab.dl.framework.mappers.FeatureMapper;
import org.campagnelab.dl.framework.mappers.MappedDimensions;
import org.nd4j.linalg.api.ndarray.INDArray;

import java.util.ArrayList;
import java.util.function.Function;

/**
 * Feature mapper that takes in a feature mapper that maps two-dimensional features, and removes masked elements in
 * between non-masked elements.
 * For example, if a feature mapper maps the following sequences, with '-' representing elements that would be masked:
 * - ACTGACTAC
 * - ACT-AC-A-
 * - AC--A----
 * The TwoDimensionalRemoveMaskFeatureMapper will map these equivalently to
 * - ACTGACTAC
 * - ACTACA---
 * - ACA------
 * Created by joshuacohen on 12/12/16.
 */
public class TwoDimensionalRemoveMaskFeatureMapper<RecordType> implements FeatureMapper<RecordType> {
    private FeatureMapper<RecordType> delegate;
    private MappedDimensions dim;
    private int featuresPerTimeStep;
    private int numTimeSteps;
    private Function<RecordType, Integer> recordToPaddingLength;


    private int[] mapperIndices = new int[]{0, 0, 0};
    private int[] maskerIndices = new int[]{0, 0};
    private ArrayList<Bounds> boundsList;

    /**
     * Creates a new feature mapper from an existing 2D feature mapper, where masked elements are removed
     * starting from position 0
     * @param delegate Delegate two-dimensional feature mapper
     */
    public TwoDimensionalRemoveMaskFeatureMapper(FeatureMapper<RecordType> delegate) {
        this(delegate, 0);
    }

    /**
     * Creates a new feature mapper from an existing 2D feature mapper, where masked elements are removed
     * starting from position startIndex
     * @param delegate Delegate two-dimensional feature mapper
     * @param startIndex starting position of where to start removing masked elements
     */
    public TwoDimensionalRemoveMaskFeatureMapper(FeatureMapper<RecordType> delegate, int startIndex) {
        this(delegate, record -> startIndex);
    }

    /**
     * Creates a new feature mapper from an existing 2D feature mapper, where masked elements are removed
     * for a given record starting from the position obtained by applying recordToPaddingLength
     * on that record
     * @param delegate Delegate two-dimensional feature mapper
     * @param recordToPaddingLength Function that, for each record, returns starting position of where to
     *                              start removing masked elements
     */
    public TwoDimensionalRemoveMaskFeatureMapper(FeatureMapper<RecordType> delegate,
                                                 Function<RecordType, Integer> recordToPaddingLength) {
        dim = delegate.dimensions();
        if (dim.numDimensions() != 2) {
            throw new RuntimeException("Delegate mapper must be two dimensional");
        }
        featuresPerTimeStep = dim.numElements(1);
        numTimeSteps = dim.numElements(2);
        this.delegate = delegate;
        this.recordToPaddingLength = recordToPaddingLength;
    }



    @Override
    public int numberOfFeatures() {
        return delegate.numberOfFeatures();
    }

    @Override
    public MappedDimensions dimensions() {
        return delegate.dimensions();
    }

    @Override
    public void prepareToNormalize(RecordType record, int indexOfRecord) {
        delegate.prepareToNormalize(record, indexOfRecord);
        boundsList = new ArrayList<>();
        if (delegate.hasMask()) {
            boolean prevMasked = true;
            Bounds currBounds = new Bounds();
            int prevBoundsIndex = -1;
            for (int i = recordToPaddingLength.apply(record); i < dim.numElements(2); i++) {
                int featureIndex = i * featuresPerTimeStep;
                if (delegate.isMasked(record, featureIndex)) {
                    if (!prevMasked) {
                        currBounds.setEnd(i);
                        if (prevBoundsIndex >= 0) {
                            currBounds.setShiftedSize(boundsList.get(prevBoundsIndex));
                        }
                        boundsList.add(currBounds);
                        prevBoundsIndex++;
                        currBounds = new Bounds();
                    }
                    prevMasked = true;
                } else {
                    if (prevMasked) {
                        currBounds.setStart(i);
                    }
                    prevMasked = false;
                }
            }
        }
    }

    @Override
    public void mapFeatures(RecordType record, INDArray inputs, int indexOfRecord) {
        mapperIndices[0] = indexOfRecord;
        for (int i = 0; i < numTimeSteps; i++) {
            for (int j = 0; j < featuresPerTimeStep; j++) {
                int featureIndex = i * featuresPerTimeStep + j;
                mapperIndices[1] = j;
                mapperIndices[2] = i;
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
        maskerIndices[0] = indexOfRecord;
        for (int i = 0; i < numTimeSteps; i++) {
            int featureIndex = i * featuresPerTimeStep;
            maskerIndices[1] = i;
            mask.putScalar(maskerIndices, isMasked(record, featureIndex) ? 1F : 0F);
        }
    }

    @Override
    public boolean isMasked(RecordType record, int featureIndex) {
        int timeStepIndex = featureIndex / featuresPerTimeStep;
        int timeStepFeatureIndex = featureIndex % featuresPerTimeStep;
        int shiftSize = 0;
        for (Bounds bounds : boundsList) {
            if (bounds.contains(timeStepIndex)) {
                shiftSize += bounds.size();
            }
        }
        int newTimeStepIndex = timeStepIndex + shiftSize;
        if (newTimeStepIndex < numTimeSteps) {
            int newFeatureIndex = newTimeStepIndex * featuresPerTimeStep + timeStepFeatureIndex;
            return delegate.isMasked(record, newFeatureIndex);
        }
        return false;
    }

    @Override
    public float produceFeature(RecordType record, int featureIndex) {
        int timeStepIndex = featureIndex / featuresPerTimeStep;
        int timeStepFeatureIndex = featureIndex % featuresPerTimeStep;
        int shiftSize = 0;
        for (Bounds bounds : boundsList) {
            if (bounds.contains(timeStepIndex)) {
                shiftSize += bounds.size();
            }
        }
        int newTimeStepIndex = timeStepIndex + shiftSize;
        if (newTimeStepIndex < numTimeSteps) {
            int newFeatureIndex = newTimeStepIndex * featuresPerTimeStep + timeStepFeatureIndex;
            return delegate.produceFeature(record, newFeatureIndex);
        }
        return 0F;
    }
}
