package org.campagnelab.dl.genotype.learning.domains.predictions;

import org.campagnelab.dl.framework.domains.prediction.PredictionInterpreter;
import org.campagnelab.dl.varanalysis.protobuf.BaseInformationRecords;
import org.nd4j.linalg.api.ndarray.INDArray;

/**
 * Created by rct66 on 11/12/16.
 */
public class HomozygousInterpreter implements PredictionInterpreter<BaseInformationRecords.BaseInformation, HomozygousPrediction> {

    double maxProbability;

    public HomozygousInterpreter(){}


    @Override
    public HomozygousPrediction interpret(BaseInformationRecords.BaseInformation record, INDArray output) {
        HomozygousPrediction pred = new HomozygousPrediction();
        pred.inspectRecord(record);
        pred.predictedHomozygousGenotype = getHomozygousPrediction(record,output);
        pred.probability = maxProbability;
        return pred;
    }

    public String getHomozygousPrediction(BaseInformationRecords.BaseInformation currentRecord, INDArray output){
        StringBuffer genotype = new StringBuffer();
        maxProbability = -1;
        int maxIndex = -1;
        for (int i = 0; i < output.length(); i++) {
            if (maxProbability < output.getDouble(0,i)){
                maxIndex = i;
                maxProbability = output.getDouble(0,i);
            }
        }
        try {
            return currentRecord.getSamples(0).getCounts(maxIndex).getToSequence();
        } catch (NullPointerException e){
            //handle non-homozygous case
            return "";
        }
    }
}
