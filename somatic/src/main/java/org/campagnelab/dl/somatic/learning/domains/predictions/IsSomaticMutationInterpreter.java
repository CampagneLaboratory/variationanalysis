package org.campagnelab.dl.somatic.learning.domains.predictions;

import org.campagnelab.dl.framework.domains.prediction.PredictionInterpreter;
import org.campagnelab.dl.varanalysis.protobuf.BaseInformationRecords;
import org.nd4j.linalg.api.ndarray.INDArray;

/**
 * Created by fac2003 on 11/12/16.
 */
public class IsSomaticMutationInterpreter implements PredictionInterpreter<BaseInformationRecords.BaseInformation, IsMutatedPrediction> {


    @Override
    public IsMutatedPrediction interpret(INDArray trueLabels, INDArray output, int exampleIndex) {
        IsMutatedPrediction prediction=new IsMutatedPrediction();
        prediction.trueLabelYes = trueLabels.getDouble(exampleIndex, 0);
        prediction.predictedLabelNo = output.getDouble(exampleIndex, 1);
        prediction.predictedLabelYes = output.getDouble(exampleIndex, 0);
        return  prediction;
    }

    @Override
    public IsMutatedPrediction interpret(BaseInformationRecords.BaseInformation record, INDArray output) {
        IsMutatedPrediction pred = new IsMutatedPrediction();
        pred.inspectRecord(record);
        pred. predictedLabelNo= (float) output.getDouble(0, 1);
        pred.predictedLabelYes = (float) output.getDouble(0, 0);
        return pred;
    }
}
