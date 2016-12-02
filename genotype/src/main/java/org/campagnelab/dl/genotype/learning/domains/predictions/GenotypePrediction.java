package org.campagnelab.dl.genotype.learning.domains.predictions;

import org.campagnelab.dl.framework.domains.prediction.BinaryClassPrediction;
import org.campagnelab.dl.varanalysis.protobuf.BaseInformationRecords;

/**
 * Created by fac2003 on 11/12/16.
 */
public class GenotypePrediction extends BinaryClassPrediction {

    int genotypeIndex;

    GenotypePrediction(int genotypeIndex){
        this.genotypeIndex = genotypeIndex;
    }

    public <BaseInformation> void inspectRecord(BaseInformationRecords.BaseInformation currentRecord) {

        trueLabelYes = (currentRecord.getSamples(0).getCountsCount() > genotypeIndex) ? (currentRecord.getSamples(0).getCounts(genotypeIndex).getIsCalled()?1d:0d) : 0d;

    }
}
