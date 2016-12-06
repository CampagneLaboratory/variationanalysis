package org.campagnelab.dl.genotype.learning.domains.predictions;

import org.campagnelab.dl.framework.domains.prediction.PredictionInterpreter;
import org.campagnelab.dl.varanalysis.protobuf.BaseInformationRecords;
import org.nd4j.linalg.api.ndarray.INDArray;

/**
 * Created by rct66 on 11/12/16.
 */
public class HomozygousInterpreter implements PredictionInterpreter<BaseInformationRecords.BaseInformation, HomozygousPrediction> {


    public HomozygousInterpreter(){}

    @Override
    public HomozygousPrediction interpret(BaseInformationRecords.BaseInformation record, INDArray output) {
        HomozygousPrediction pred = new HomozygousPrediction();
        pred.inspectRecord(record);
        pred.predictedProbabilities = new double[output.length()];
        for (int i = 0; i < pred.predictedProbabilities.length; i++){
            pred.predictedProbabilities[i] = output.getDouble(0,i);
        }
        return pred;
    }
}
