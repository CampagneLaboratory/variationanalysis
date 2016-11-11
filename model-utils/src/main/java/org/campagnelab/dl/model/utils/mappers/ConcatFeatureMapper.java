package org.campagnelab.dl.model.utils.mappers;

import org.campagnelab.dl.varanalysis.protobuf.BaseInformationRecords;
import org.nd4j.linalg.api.ndarray.INDArray;

import java.util.Arrays;

/**
 * Concatenate features from different mappers.
 * Created by fac2003 on 5/24/16.
 */
public class ConcatFeatureMapper<RecordType> implements FeatureMapper<RecordType> {

    protected FeatureMapper<RecordType>[] mappers;
    protected int numFeatures = 0;
    protected int[] offsets;
    private boolean normalizedCalled;

    @SafeVarargs
    public ConcatFeatureMapper(FeatureMapper<RecordType> ... featureMappers) {

        this.mappers = featureMappers;
        int offset = 0;
        int i = 1;
        offsets = new int[featureMappers.length + 1];
        offsets[0] = 0;
        for (FeatureMapper<RecordType> calculator : mappers) {
            numFeatures += calculator.numberOfFeatures();
            ;
            offsets[i] = numFeatures;

            i++;
        }
    }

    @Override
    public int numberOfFeatures() {

        return numFeatures;
    }

    @Override
    public void prepareToNormalize(RecordType record, int indexOfRecord) {
        for (FeatureMapper<RecordType>  calculator : mappers) {
            calculator.prepareToNormalize(record, indexOfRecord);
        }
        normalizedCalled=true;
    }


    @Override
    public void mapFeatures(RecordType record, INDArray inputs, int indexOfRecord) {
        assert normalizedCalled :"prepareToNormalize must be called before mapFeatures.";
        int offset = 0;
        final int[] indicesOuter = {0, 0};
        for (FeatureMapper<RecordType>  delegate : mappers) {

            final int delNumFeatures = delegate.numberOfFeatures();
            for (int j = 0; j < delNumFeatures; j++) {
                indicesOuter[0] = indexOfRecord;
                indicesOuter[1] = j + offset;
                inputs.putScalar(indicesOuter, delegate.produceFeature(record, j));
            }
            offset += delNumFeatures;
        }
    }

    @Override
    public float produceFeature(RecordType record, int featureIndex) {
        int indexOfDelegate = Arrays.binarySearch(offsets, featureIndex);
        if (indexOfDelegate < 0) {
            indexOfDelegate = -(indexOfDelegate + 1) - 1;
        }
        return this.mappers[indexOfDelegate].produceFeature(record, featureIndex - offsets[indexOfDelegate]);
    }



}
