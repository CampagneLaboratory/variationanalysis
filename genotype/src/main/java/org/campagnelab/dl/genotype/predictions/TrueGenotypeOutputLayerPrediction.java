package org.campagnelab.dl.genotype.predictions;

import org.campagnelab.dl.framework.domains.prediction.Prediction;
import org.campagnelab.dl.varanalysis.protobuf.BaseInformationRecords;

/**
 * Created by joshuacohen on 2/21/17.
 */
public class TrueGenotypeOutputLayerPrediction extends Prediction {
    /**
     * Actual genotype
     */
    public String trueGenotype;
    /**
     * Genotype called by the model
     */
    public String predictedGenotype;
    /**
     * The probability of the called genotype according to the model. Forecast probability.
     */
    public double overallProbability;
    /**
     * If the predicted genotype corresponds to an indel
     */
    public boolean isPredictedIndel;

    public void inspectRecord(BaseInformationRecords.BaseInformation record) {
        trueGenotype = record.getTrueGenotype();
    }
}
