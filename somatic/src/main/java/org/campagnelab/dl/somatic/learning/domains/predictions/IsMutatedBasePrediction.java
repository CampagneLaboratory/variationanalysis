package org.campagnelab.dl.somatic.learning.domains.predictions;

import org.campagnelab.dl.varanalysis.protobuf.BaseInformationRecords;

/**
 * Created by fac2003 on 12/5/16.
 */
public class IsMutatedBasePrediction extends IsMutatedPrediction {
    public <BaseInformation> void inspectRecord(BaseInformationRecords.BaseInformation currentRecord) {
        super.inspectRecord(currentRecord);

    }
}
