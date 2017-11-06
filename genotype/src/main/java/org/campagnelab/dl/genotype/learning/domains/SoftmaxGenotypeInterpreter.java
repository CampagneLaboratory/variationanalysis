package org.campagnelab.dl.genotype.learning.domains;

import org.campagnelab.dl.framework.domains.prediction.PredictionInterpreter;
import org.campagnelab.dl.genotype.helpers.GenotypeHelper;
import org.campagnelab.dl.genotype.predictions.SoftmaxGenotypePrediction;
import org.campagnelab.dl.varanalysis.protobuf.BaseInformationRecords;
import org.nd4j.linalg.api.ndarray.INDArray;

/**
 * Interprets the output of SoftmaxGenotype layer.
 * Created by fac2003 on 2/7/17.
 */
public class SoftmaxGenotypeInterpreter implements PredictionInterpreter<BaseInformationRecords.BaseInformation,
        SoftmaxGenotypePrediction> {

    int maxCalledAlleles;
    private double maxProbability;
    private int numBits;

    /**
     * Construct the interpreter.
     * @param maxCalledAlleles The maximum number of possible called alleles.
     */
    public SoftmaxGenotypeInterpreter(int maxCalledAlleles) {
        this.maxCalledAlleles = maxCalledAlleles;
        numBits = (int) Math.pow(2, maxCalledAlleles);
    }

    @Override
    public SoftmaxGenotypePrediction interpret(INDArray trueLabels, INDArray output, int predictionIndex) {
        SoftmaxGenotypePrediction result = new SoftmaxGenotypePrediction();

        result.predictedGenotypeIndex = readPredicted(output, result, predictionIndex);
        result.probability = maxProbability;
        result.trueGenotypeIndex = readPredicted(trueLabels, result, predictionIndex);
        result.numBits=numBits;
        return result;
    }

    @Override
    public SoftmaxGenotypePrediction interpret(BaseInformationRecords.BaseInformation record, INDArray output) {
        SoftmaxGenotypePrediction result = new SoftmaxGenotypePrediction();
        result.numBits=numBits;
        final String trueGenotype = record.getTrueGenotype();
        result.trueGenotype = trueGenotype;

        // interpret the prediction
        result.predictedGenotypeIndex = readPredicted(output, result, 0);
        result.probability = maxProbability;
        return result;
    }

    private int readPredicted(INDArray output, SoftmaxGenotypePrediction result, int predictionIndex) {
        maxProbability = -1;
        int maxIndex = -1;

        for (int i = 0; i < numBits; i++) {
            double outputDouble = output.getDouble(predictionIndex, i);
            if (maxProbability < outputDouble) {
                maxIndex = i;
                maxProbability = outputDouble;
            }
        }
        return maxIndex;
    }
}
