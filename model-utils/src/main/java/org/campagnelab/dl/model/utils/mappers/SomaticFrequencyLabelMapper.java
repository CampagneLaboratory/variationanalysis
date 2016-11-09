package org.campagnelab.dl.model.utils.mappers;

import org.campagnelab.dl.varanalysis.protobuf.BaseInformationRecords;
import org.nd4j.linalg.api.ndarray.INDArray;

/**
 * Label: frequency of somatic mutation.
 * Created by fac2003 on 11/8/16.
 */
public class SomaticFrequencyLabelMapper implements LabelMapper<BaseInformationRecords.BaseInformation> {
    @Override
    public int numberOfLabels() {
        return 1;
    }

    @Override
    public void mapLabels(BaseInformationRecords.BaseInformation record, INDArray labels, int indexOfRecord) {
        labels.putScalar(indexOfRecord, record.getFrequencyOfMutation());
    }

    @Override
    public float produceLabel(BaseInformationRecords.BaseInformation record, int labelIndex) {
        assert labelIndex == 0 || labelIndex == 1 : "only one label.";

        if (labelIndex == 0)
            return record.getFrequencyOfMutation();
        else {
            return 1 - record.getFrequencyOfMutation();
        }

    }
}
