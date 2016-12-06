package org.campagnelab.dl.genotype.learning.domains.predictions;

import org.campagnelab.dl.framework.domains.prediction.PredictionInterpreter;
import org.campagnelab.dl.somatic.learning.domains.predictions.IsMutatedPrediction;
import org.campagnelab.dl.varanalysis.protobuf.BaseInformationRecords;
import org.nd4j.linalg.api.ndarray.INDArray;

/**
 * Created by rct66 on 11/12/16.
 */
public class GenotypeInterpreter implements PredictionInterpreter<BaseInformationRecords.BaseInformation, GenotypePrediction> {

    int genotypeIndex;

    public GenotypeInterpreter(int genotypeIndex){
        this.genotypeIndex = genotypeIndex;
    }

    @Override
    public GenotypePrediction interpret(BaseInformationRecords.BaseInformation record, INDArray output) {
        GenotypePrediction pred = new GenotypePrediction(genotypeIndex);
        pred.inspectRecord(record);
        pred.predictedLabelYes = (float) output.getDouble(0, 0);
        pred.predictedLabelNo= (float) output.getDouble(0, 1);
        return pred;
    }
}
