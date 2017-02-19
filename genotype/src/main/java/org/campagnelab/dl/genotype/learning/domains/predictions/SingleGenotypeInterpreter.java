package org.campagnelab.dl.genotype.learning.domains.predictions;

import org.campagnelab.dl.framework.domains.prediction.PredictionInterpreter;
import org.campagnelab.dl.varanalysis.protobuf.BaseInformationRecords;
import org.nd4j.linalg.api.ndarray.INDArray;

/**
 * Interpret a predicton for a single genotype.
 * Created by rct66 on 11/12/16.
 */
public class SingleGenotypeInterpreter extends SortingCountInterpreter<SingleGenotypePrediction>
        implements PredictionInterpreter<BaseInformationRecords.BaseInformation, SingleGenotypePrediction> {

    int genotypeIndex;

    public SingleGenotypeInterpreter(int genotypeIndex, boolean sort) {
        super(sort);
        this.genotypeIndex = genotypeIndex;
    }

    ;


    @Override
    public SingleGenotypePrediction interpret(INDArray trueLabels, INDArray output, int predictionIndex) {
        SingleGenotypePrediction pred = new SingleGenotypePrediction();
        try {
            pred.predictedSingleGenotype = Integer.toString(genotypeIndex);
            if (genotypeIndex > 4) {
                pred.isPredictedIndel = true;
            }
        } catch (IndexOutOfBoundsException e) {
            pred.predictedSingleGenotype = ".";
        }
        pred.probabilityIsCalled = output.getDouble(predictionIndex, 0);
        pred.trueIsCalled = trueLabels.getDouble(predictionIndex, 0) > 0.5;
        return pred;
    }

    @Override
    public SingleGenotypePrediction interpret(BaseInformationRecords.BaseInformation record, INDArray output) {
        SingleGenotypePrediction pred = new SingleGenotypePrediction();
        int predictionIndex = 0;

        try {
            BaseInformationRecords.CountInfo counts = sort(record).getSamples(0).getCounts(indexPermutation[genotypeIndex]);
            pred.predictedSingleGenotype = counts.getToSequence();
            pred.isPredictedIndel = counts.getIsIndel();
        } catch (IndexOutOfBoundsException e) {
            pred.predictedSingleGenotype = ".";
        }
        pred.probabilityIsCalled = output.getDouble(predictionIndex, 0);
        return pred;
    }


}
