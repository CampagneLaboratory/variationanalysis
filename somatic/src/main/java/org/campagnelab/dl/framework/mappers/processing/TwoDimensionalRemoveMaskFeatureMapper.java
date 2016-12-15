package org.campagnelab.dl.framework.mappers.processing;

import org.campagnelab.dl.framework.mappers.FeatureMapper;
import org.campagnelab.dl.framework.mappers.MappedDimensions;
import org.nd4j.linalg.api.ndarray.INDArray;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;

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
    static private Logger LOG = LoggerFactory.getLogger(TwoDimensionalRemoveMaskFeatureMapper.class);

    private FeatureMapper<RecordType> delegate;
    private MappedDimensions dim;
    private HashSet<RecordType> normalizedCalled;
    private HashMap<RecordType, ArrayList<Bounds>> boundsMap;
    private int featuresPerTimeStep;
    private int numTimeSteps;
    private int startIndex;


    private int[] mapperIndices = new int[]{0, 0, 0};
    private int[] maskerIndices = new int[]{0, 0};


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
        dim = delegate.dimensions();
        if (dim.numDimensions() != 2) {
            throw new RuntimeException("Delegate mapper must be two dimensional");
        }
        featuresPerTimeStep = dim.numElements(1);
        numTimeSteps = dim.numElements(2);
        this.delegate = delegate;
        this.startIndex = startIndex;
        boundsMap = new HashMap<>();
        normalizedCalled = new HashSet<>();
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
        if (indexOfRecord == -1) {
            LOG.warn("prepareToNormalize should be called before mapping/masking or calling produceFeature/isMasked");
        }
        if (delegate.hasMask()) {
            boolean prevMasked = true;
            Bounds currBounds = new Bounds();
            ArrayList<Bounds> boundsList = new ArrayList<>();
            int prevBoundsIndex = -1;
            for (int i = startIndex; i < dim.numElements(2); i++) {
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
            boundsMap.put(record, boundsList);
        }
        normalizedCalled.add(record);
    }

    @Override
    public void mapFeatures(RecordType record, INDArray inputs, int indexOfRecord) {
        mapperIndices[0] = indexOfRecord;
        for (int i = 0; i < dim.numElements(2); i++) {
            for (int j = 0; j < dim.numElements(1); j++) {
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
        for (int i = 0; i < dim.numElements(2); i++) {
            int featureIndex = i * featuresPerTimeStep;
            maskerIndices[1] = i;
            mask.putScalar(maskerIndices, isMasked(record, featureIndex) ? 1F : 0F);
        }
    }

    @Override
    public boolean isMasked(RecordType record, int featureIndex) {
        if (!normalizedCalled.contains(record)) {
            prepareToNormalize(record, -1);
        }
        int timeStepIndex = featureIndex / featuresPerTimeStep;
        int timeStepFeatureIndex = featureIndex % featuresPerTimeStep;
        ArrayList<Bounds> boundsList = boundsMap.get(record);
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
        if (!normalizedCalled.contains(record)) {
            prepareToNormalize(record, -1);
        }
        int timeStepIndex = featureIndex / featuresPerTimeStep;
        int timeStepFeatureIndex = featureIndex % featuresPerTimeStep;
        ArrayList<Bounds> boundsList = boundsMap.get(record);
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
