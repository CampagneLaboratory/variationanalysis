package org.campagnelab.dl.genotype.learning.domains.predictions;

import org.campagnelab.dl.framework.domains.prediction.PredictionInterpreter;
import org.campagnelab.dl.varanalysis.protobuf.BaseInformationRecords;
import org.nd4j.linalg.api.ndarray.INDArray;

/**
 * Created by rct66 on 11/12/16.
 */
public class SingleGenotypeInterpreter implements PredictionInterpreter<BaseInformationRecords.BaseInformation, SingleGenotypePrediction> {

    int genotypeIndex;

    public SingleGenotypeInterpreter(int genotypeIndex) {
        this.genotypeIndex = genotypeIndex;

    };


    @Override
    public SingleGenotypePrediction interpret(BaseInformationRecords.BaseInformation record, INDArray output) {
        SingleGenotypePrediction pred = new SingleGenotypePrediction();
        pred.predictedSingleGenotype = record.getSamples(0).getCounts(genotypeIndex).getToSequence();
        pred.probability = output.getDouble(0,0);
        return pred;
    }


}
