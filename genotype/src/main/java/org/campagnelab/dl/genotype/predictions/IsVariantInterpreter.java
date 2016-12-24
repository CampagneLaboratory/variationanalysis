package org.campagnelab.dl.genotype.predictions;

import org.campagnelab.dl.framework.domains.prediction.PredictionInterpreter;
import org.campagnelab.dl.framework.mappers.BooleanLabelMapper;
import org.campagnelab.dl.varanalysis.protobuf.BaseInformationRecords;
import org.nd4j.linalg.api.ndarray.INDArray;

/**
 * This interpreter extracts isVariant predictions (confidence in a genotype being a true variant).
 * Created by fac2003 on 12/22/16.
 */
public class IsVariantInterpreter implements PredictionInterpreter<BaseInformationRecords.BaseInformation,
        IsVariantPrediction> {

    @Override
    public IsVariantPrediction interpret(INDArray trueLabels, INDArray output, int predictionIndex) {
        IsVariantPrediction p = new IsVariantPrediction();
        p.isVariantPredicted = output.getDouble(predictionIndex, BooleanLabelMapper.IS_TRUE) > 0.5;
        p.isVariantTruth = trueLabels.getDouble(predictionIndex, BooleanLabelMapper.IS_TRUE) > 0.5;
        p.probability =Math.max(
                output.getDouble(predictionIndex, BooleanLabelMapper.IS_TRUE),
                output.getDouble(predictionIndex, BooleanLabelMapper.IS_FALSE))
        ;

        return p;
    }

    @Override
    public IsVariantPrediction interpret(BaseInformationRecords.BaseInformation record, INDArray output) {
        IsVariantPrediction p = new IsVariantPrediction();
        p.isVariantTruth = record.getSamples(0).getIsVariant();
        int predictionIndex = 0;
        p.isVariantPredicted = output.getDouble(predictionIndex, BooleanLabelMapper.IS_TRUE) > 0.5;

        return p;
    }
}
