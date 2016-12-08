package org.campagnelab.dl.genotype.learning.domains.predictions;

import org.campagnelab.dl.framework.domains.prediction.Prediction;
import org.campagnelab.dl.varanalysis.protobuf.BaseInformationRecords;

/**
 * Created by rct66 on 11/12/16.
 */
public class HomozygousPrediction extends Prediction {


    public String trueGenotype;
    public String predictedHomozygousGenotype;
    public double probability;

    public <BaseInformation> void inspectRecord(BaseInformationRecords.BaseInformation currentRecord) {
        trueGenotype = getGenotype(currentRecord);
        return;
    }

    public String getGenotype(BaseInformationRecords.BaseInformation currentRecord){
        StringBuffer genotype = new StringBuffer();
        for (int i = 0; i < currentRecord.getSamples(0).getCountsCount(); i++) {
            BaseInformationRecords.CountInfo count = currentRecord.getSamples(0).getCounts(i);
            if (count != null) {
                if (count.getIsCalled()){
                    genotype.append(count.getToSequence() + ",");
                }
            }
        }
        return genotype.substring(0,genotype.length());
    }
}

