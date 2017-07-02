package org.campagnelab.dl.genotype.learning.domains;

import org.campagnelab.dl.framework.domains.prediction.PredictionInterpreter;

import org.campagnelab.dl.genotype.helpers.GenotypeHelper;
import org.campagnelab.dl.genotype.predictions.GenotypePrediction;
import org.campagnelab.dl.varanalysis.protobuf.BaseInformationRecords;
import org.nd4j.linalg.api.ndarray.INDArray;

/**
 * Interprets the output of NumDistinctAllelesLabelMapper.
 * Created by fac2003 on 12/20/16.
 */
public class NumDistinctAllelesInterpreter implements PredictionInterpreter<BaseInformationRecords.BaseInformation,
        NumDistinctAllelesOutputLayerPrediction> {

    int ploidy;
    private double maxProbability;

    public NumDistinctAllelesInterpreter(int ploidy) {
        this.ploidy = ploidy;
    }

    @Override
    public NumDistinctAllelesOutputLayerPrediction interpret(INDArray trueLabels, INDArray output, int predictionIndex) {
        NumDistinctAllelesOutputLayerPrediction result = new NumDistinctAllelesOutputLayerPrediction();

        result.predictedValue = readPredicted(output, result,predictionIndex);
        result.trueValue = readPredicted(trueLabels, result,predictionIndex);
        result.probability=maxProbability;
        return result;
    }

    @Override
    public NumDistinctAllelesOutputLayerPrediction interpret(BaseInformationRecords.BaseInformation record, INDArray output) {
        NumDistinctAllelesOutputLayerPrediction result = new NumDistinctAllelesOutputLayerPrediction();
        final String trueGenotype = record.getTrueGenotype();
        result.trueValue = GenotypeHelper.getAlleles(trueGenotype).size();
        // interpret the prediction
        result.predictedValue = readPredicted(output, result,0);
        result.probability=maxProbability;
        return result;
    }

    private int readPredicted(INDArray output, NumDistinctAllelesOutputLayerPrediction result, int  predictionIndex) {
        maxProbability = -1;
        int maxIndex = -1;

        for (int i = 0; i < ploidy+1; i++) {
            double outputDouble = output.getDouble(predictionIndex, i);
            if (maxProbability < outputDouble) {
                maxIndex = i;
                maxProbability = outputDouble;
            }
        }
       return maxIndex;
    }
}
