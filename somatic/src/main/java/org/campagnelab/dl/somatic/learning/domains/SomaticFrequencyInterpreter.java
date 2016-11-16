package org.campagnelab.dl.somatic.learning.domains;

import org.campagnelab.dl.somatic.learning.domains.predictions.SomaticFrequencyPrediction;
import org.campagnelab.dl.framework.domains.prediction.PredictionInterpreter;
import org.campagnelab.dl.varanalysis.protobuf.BaseInformationRecords;
import org.nd4j.linalg.api.ndarray.INDArray;

/**
 * Created by fac2003 on 11/12/16.
 */
public class SomaticFrequencyInterpreter implements PredictionInterpreter<BaseInformationRecords.BaseInformation, SomaticFrequencyPrediction> {
    @Override
    public SomaticFrequencyPrediction interpret(BaseInformationRecords.BaseInformation record, INDArray output) {
        SomaticFrequencyPrediction pred = new SomaticFrequencyPrediction();
        // set trueValue to 0 if the record is not mutated. We need this to calculate rmse.
        pred.trueValue = record.hasFrequencyOfMutation() || record.hasMutated() ? record.getFrequencyOfMutation() : null;
        pred.predictedValue = output.getFloat(0);
        return pred;
    }
}
