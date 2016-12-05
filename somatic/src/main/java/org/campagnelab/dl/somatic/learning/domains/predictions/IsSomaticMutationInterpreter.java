package org.campagnelab.dl.somatic.learning.domains.predictions;

import org.campagnelab.dl.framework.domains.prediction.PredictionInterpreter;
import org.campagnelab.dl.varanalysis.protobuf.BaseInformationRecords;
import org.nd4j.linalg.api.ndarray.INDArray;

/**
 * Created by fac2003 on 11/12/16.
 */
public class IsSomaticMutationInterpreter implements PredictionInterpreter<BaseInformationRecords.BaseInformation, IsMutatedPrediction> {


    @Override
    public IsMutatedPrediction interpret(INDArray trueLabels, INDArray[] outputs, int predictionIndex) {
        IsMutatedPrediction prediction=new IsMutatedPrediction();
        prediction.trueLabelYes = trueLabels.getDouble(predictionIndex, 1);
        int outputIndex=0;
        prediction.predictedLabelNo = outputs[outputIndex].getDouble(predictionIndex, 0);
        prediction.predictedLabelYes = outputs[outputIndex].getDouble(predictionIndex, 1);
        return  prediction;
    }

    @Override
    public IsMutatedPrediction interpret(BaseInformationRecords.BaseInformation record, INDArray output) {
        IsMutatedPrediction pred = new IsMutatedPrediction();
        pred.inspectRecord(record);
        pred.predictedLabelYes = (float) output.getDouble(0, 0);
        pred. predictedLabelNo= (float) output.getDouble(0, 1);
        return pred;
    }
}
