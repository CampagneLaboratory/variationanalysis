package org.campagnelab.dl.varanalysis.learning.domains.predictions;

import org.campagnelab.dl.varanalysis.protobuf.BaseInformationRecords;
import org.nd4j.linalg.api.ndarray.INDArray;

/**
 * Created by fac2003 on 11/12/16.
 */
public class IsSomaticMutationInterpreter implements PredictionInterpreter<BaseInformationRecords.BaseInformation, IsMutatedPrediction> {


    @Override
    public IsMutatedPrediction interpret(BaseInformationRecords.BaseInformation record, INDArray output) {
        IsMutatedPrediction pred = new IsMutatedPrediction();
        pred.inspectRecord(record);
        pred.predictedLabelYes = (float) output.getDouble(0, 0);
        pred. predictedLabelNo= (float) output.getDouble(0, 1);
        return pred;
    }
}
