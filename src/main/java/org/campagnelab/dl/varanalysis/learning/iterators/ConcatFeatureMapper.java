package org.campagnelab.dl.varanalysis.learning.iterators;

import org.campagnelab.dl.varanalysis.format.PosRecord;
import org.nd4j.linalg.api.ndarray.INDArray;
import org.nd4j.linalg.cpu.nativecpu.NDArray;
import org.nd4j.linalg.factory.Nd4j;

import java.util.Arrays;

/**
 * Concatenate features from different mappers.
 * Created by fac2003 on 5/24/16.
 */
public class ConcatFeatureMapper implements FeatureMapper {
    private FeatureMapper mappers[];
    private int numFeatures = 0;
    private int[] offsets;

    public ConcatFeatureMapper(FeatureMapper... featureMappers) {
        this.mappers = featureMappers;
        int offset = 0;
        int i = 1;
        offsets=new int[featureMappers.length+1];
        offsets[0]=0;
        for (FeatureMapper calculator : mappers) {
            numFeatures += calculator.numberOfFeatures();
            offset += numFeatures;
            offsets[i] = offset;

            i++;
        }
    }

    @Override
    public int numberOfFeatures() {

        return numFeatures;
    }


    @Override
    public void mapFeatures(PosRecord record, INDArray inputs, int indexOfRecord) {
        int offset = 0;
        int[] indicesOuter = {0};
        for (FeatureMapper delegate : mappers) {

            final int delNumFeatures = delegate.numberOfFeatures();

            for (int j = 0; j < delNumFeatures; j++) {
                indicesOuter[0] = j + offset;
                inputs.putScalar(indicesOuter, delegate.produceFeature(record, j));
            }
            offset += delNumFeatures;
        }
    }

    @Override
    public float produceFeature(PosRecord record, int featureIndex) {
        int indexOfDelegate = Arrays.binarySearch(offsets, featureIndex);
        if (indexOfDelegate<0) {
            indexOfDelegate=-(indexOfDelegate+1)-1;
        }
        return this.mappers[indexOfDelegate].produceFeature(record, featureIndex);
    }

}
