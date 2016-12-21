package org.campagnelab.dl.somatic.learning.domains;

import org.campagnelab.dl.framework.domains.prediction.PredictionInterpreter;
import org.campagnelab.dl.somatic.learning.domains.predictions.IsMutatedPrediction;
import org.campagnelab.dl.somatic.learning.domains.predictions.IsSomaticMutationInterpreter;
import org.campagnelab.dl.varanalysis.protobuf.BaseInformationRecords;
import org.nd4j.linalg.api.ndarray.INDArray;

/**
 * Interpret whether the output indicates the site is mutated. Used with IsBaseMutatedMapper.
 * Created by fac2003 on 12/5/16.
 */
public class IsBaseMutatedInterpreter extends IsSomaticMutationInterpreter implements PredictionInterpreter<BaseInformationRecords.BaseInformation, IsMutatedPrediction> {

    @Override
    public IsMutatedPrediction interpret(INDArray trueLabels, INDArray output, int predictionIndex) {
        IsMutatedPrediction prediction = new IsMutatedPrediction();
        int outputIndex = 0; // first output is IsBaseMutated
        prediction.predictedLabelNo = output.getDouble(predictionIndex, 0);
        prediction.predictedLabelYes = 1f -output.getDouble(predictionIndex, 0);
        prediction.trueLabelYes = 1-trueLabels.getDouble(predictionIndex, 0);
        return prediction;
    }

    @Override
    public IsMutatedPrediction interpret(BaseInformationRecords.BaseInformation record, INDArray output) {
        IsMutatedPrediction prediction = new IsMutatedPrediction();
        prediction.inspectRecord(record);
        int predictionIndex=0;
        prediction.predictedLabelNo = output.getDouble(predictionIndex, 0);
        prediction.predictedLabelYes = 1f- output.getDouble(predictionIndex, 0);
        // TODO extract information about the sorted genotype called using the next 5 floats in the output.
        return prediction;
    }
}
