package org.campagnelab.dl.framework.mappers;

import it.unimi.dsi.fastutil.ints.IntArraySet;
import org.nd4j.linalg.api.ndarray.INDArray;

import java.util.Arrays;

/**
 * Concatenate labels from different mappers.
 * Created by fac2003 on 12/5/16.
 */
public class ConcatLabelMapper<RecordType> implements LabelMapper<RecordType> {

    protected LabelMapper<RecordType>[] mappers;
    protected int numFeatures = 0;
    protected int[] offsets;
    private boolean normalizedCalled;

    @SafeVarargs
    public ConcatLabelMapper(LabelMapper<RecordType>... labelMappers) {
        IntArraySet dimensions = new IntArraySet();

        this.mappers = labelMappers;
        int offset = 0;
        int i = 1;
        offsets = new int[labelMappers.length + 1];
        offsets[0] = 0;
        for (LabelMapper<RecordType> calculator : mappers) {
            numFeatures += calculator.numberOfLabels();
            offsets[i] = numFeatures;
            dimensions.add(calculator.dimensions().numDimensions());
            i++;
        }
        assert labelMappers.length==0 || dimensions.size()==1: "All feature mappers must have the same dimensions to be concatenated.";
    }

    @Override
    public int numberOfLabels() {

        return numFeatures;
    }

    @Override
    public MappedDimensions dimensions() {
        return new MappedDimensions(numberOfLabels());
    }

    @Override
    public void prepareToNormalize(RecordType record, int indexOfRecord) {
        for (LabelMapper<RecordType> calculator : mappers) {
            calculator.prepareToNormalize(record, indexOfRecord);
        }
        normalizedCalled = true;
    }


    @Override
    public void mapLabels(RecordType record, INDArray inputs, int indexOfRecord) {
        assert normalizedCalled : "prepareToNormalize must be called before mapFeatures.";
        int offset = 0;
        final int[] indicesOuter = {0, 0};
        for (LabelMapper<RecordType> delegate : mappers) {

            final int delNumFeatures = delegate.numberOfLabels();
            for (int j = 0; j < delNumFeatures; j++) {
                indicesOuter[0] = indexOfRecord;
                indicesOuter[1] = j + offset;
                inputs.putScalar(indicesOuter, delegate.produceLabel(record, j));
            }
            offset += delNumFeatures;
        }
    }

    @Override
    public boolean hasMask() {
        boolean requiresMask = false;
        for (LabelMapper<RecordType> calculator : mappers) {
            requiresMask |= calculator.hasMask();
        }
        return requiresMask;
    }

    @Override
    public void maskLabels(RecordType record, INDArray mask, int indexOfRecord) {
        if (hasMask()) {
            int offset = 0;
            final int[] indicesOuter = {0, 0};
            for (LabelMapper<RecordType> delegate : mappers) {

                final int delNumFeatures = delegate.numberOfLabels();
                for (int j = 0; j < delNumFeatures; j++) {
                    indicesOuter[0] = indexOfRecord;
                    indicesOuter[1] = j + offset;
                    mask.putScalar(indicesOuter, delegate.isMasked(record, j) ? 1 : 0);
                }
                offset += delNumFeatures;
            }
        }
    }

    @Override
    public boolean isMasked(RecordType record, int featureIndex) {
        int indexOfDelegate = Arrays.binarySearch(offsets, featureIndex);
        if (indexOfDelegate < 0) {
            indexOfDelegate = -(indexOfDelegate + 1) - 1;
        }
        return this.mappers[indexOfDelegate].isMasked(record, featureIndex - offsets[indexOfDelegate]);
    }

    @Override
    public float produceLabel(RecordType record, int featureIndex) {
        int indexOfDelegate = Arrays.binarySearch(offsets, featureIndex);
        if (indexOfDelegate < 0) {
            indexOfDelegate = -(indexOfDelegate + 1) - 1;
        }
        return this.mappers[indexOfDelegate].produceLabel(record, featureIndex - offsets[indexOfDelegate]);
    }


}
