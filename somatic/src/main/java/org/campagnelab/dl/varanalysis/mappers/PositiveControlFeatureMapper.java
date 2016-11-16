package org.campagnelab.dl.varanalysis.mappers;

import org.campagnelab.dl.varanalysis.protobuf.BaseInformationRecords;
import org.nd4j.linalg.api.ndarray.INDArray;

/**
 * The FeatureMapper to test for the second iteration.
 * Created by fac2003 on 5/24/16.
 */
public class PositiveControlFeatureMapper extends NoMaskFeatureMapper<BaseInformationRecords.BaseInformationOrBuilder> {

    @Override
    public int numberOfFeatures() {
        return 2;
    }

    @Override
    public void prepareToNormalize(BaseInformationRecords.BaseInformationOrBuilder record, int indexOfRecord) {

    }

    @Override
    public void mapFeatures(BaseInformationRecords.BaseInformationOrBuilder record, INDArray inputs, int indexOfRecord) {
       for (int i=0;i<numberOfFeatures();i++) {
           inputs.putScalar(new int[]{indexOfRecord, i}, produceFeature(record, 0));
           //       System.out.println(inputs);
       }
    }

    @Override
    public float produceFeature(BaseInformationRecords.BaseInformationOrBuilder record, int featureIndex) {
        return record.getMutated() ? 1 : 0;
    }

}
