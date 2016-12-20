package org.campagnelab.dl.genotype.learning.domains.predictions;

import org.campagnelab.dl.framework.domains.prediction.PredictionInterpreter;
import org.campagnelab.dl.varanalysis.protobuf.BaseInformationRecords;
import org.nd4j.linalg.api.ndarray.INDArray;

/**
 * Created by rct66 on 11/12/16.
 */
public class SingleGenotypeInterpreter extends  SortingCountInterpreter<SingleGenotypePrediction>
        implements PredictionInterpreter<BaseInformationRecords.BaseInformation, SingleGenotypePrediction> {

    int genotypeIndex;

    public SingleGenotypeInterpreter(int genotypeIndex, boolean sort) {
        super(sort);
        this.genotypeIndex = genotypeIndex;

    };


    @Override
    public SingleGenotypePrediction interpret(INDArray trueLabels, INDArray[] outputs, int predictionIndex) {
        throw new RuntimeException("a wrong interpret method was called on the homozygous interpeter");
    }

    @Override
    public SingleGenotypePrediction interpret(BaseInformationRecords.BaseInformation record, INDArray output) {
        SingleGenotypePrediction pred = new SingleGenotypePrediction();
        try {
            pred.predictedSingleGenotype = sort(record).getSamples(0).getCounts(genotypeIndex).getToSequence();
        } catch (IndexOutOfBoundsException e) {
            pred.predictedSingleGenotype = ".";
        }
        pred.probabilityIsCalled = output.getDouble(0,0);
        return pred;
    }


}
