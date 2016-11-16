package org.campagnelab.dl.varanalysis.mappers;

import org.campagnelab.dl.varanalysis.protobuf.BaseInformationRecords;
import org.nd4j.linalg.api.ndarray.INDArray;

/**
 * Label: frequency of somatic mutation.
 * Created by fac2003 on 11/8/16.
 */
public class SomaticFrequencyLabelMapper extends NoMasksLabelMapper<BaseInformationRecords.BaseInformation> {
    @Override
    public int numberOfLabels() {
        return 1;
    }

    int[] indices = new int[]{0, 0};

    @Override
    public void mapLabels(BaseInformationRecords.BaseInformation record, INDArray labels, int indexOfRecord) {
        indices[0] = indexOfRecord;

        for (int labelIndex = 0; labelIndex < numberOfLabels(); labelIndex++) {
            indices[1] = labelIndex;
            labels.putScalar(indices, produceLabel(record, labelIndex));
        }
    }

    @Override
    public float produceLabel(BaseInformationRecords.BaseInformation record, int labelIndex) {
        assert labelIndex == 0 || labelIndex == 1 : "only one label.";

        return record.getFrequencyOfMutation();
    }


}
