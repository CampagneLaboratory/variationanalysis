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

    int sortedGenotypeIndex;

    public SingleGenotypeInterpreter(int sortedGenotypeIndex, boolean sort) {
        super(sort);
        this.sortedGenotypeIndex = sortedGenotypeIndex;
    }

    @Override
    public SingleGenotypePrediction interpret(INDArray trueLabels, INDArray output, int exampleIndex) {
        SingleGenotypePrediction pred = new SingleGenotypePrediction();
        pred.sortedCountIndex = sortedGenotypeIndex;
        pred.probabilityIsCalled = output.getDouble(exampleIndex, 0);
        pred.trueIsCalled = trueLabels.getDouble(exampleIndex, 0) > 0.5;
        return pred;
    }

    @Override
    public SingleGenotypePrediction interpret(BaseInformationRecords.BaseInformation record, INDArray output) {
        SingleGenotypePrediction pred = new SingleGenotypePrediction();
        int predictionIndex = 0;
        pred.sortedCountIndex = sortedGenotypeIndex;
        pred.probabilityIsCalled = output.getDouble(predictionIndex, 0);
        sortedCountRecord = sort(record);
        if (sortedGenotypeIndex < sortedCountRecord.getSamples(0).getCountsCount()) {
            // ONLY so many genotypes stored.
            pred.trueIsCalled = sortedCountRecord.getSamples(0).getCounts(sortedGenotypeIndex).getIsCalled();
        } else {
            pred.trueIsCalled = false;
        }
        return pred;
    }
}
