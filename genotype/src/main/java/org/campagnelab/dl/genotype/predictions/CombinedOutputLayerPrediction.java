package org.campagnelab.dl.genotype.predictions;

import org.campagnelab.dl.framework.domains.prediction.Prediction;
import org.campagnelab.dl.varanalysis.protobuf.BaseInformationRecords;

/**
 * A class to store information interpreted from the Combined layer output.
 * Created by fac2003 on 12/22/16.
 */
public class CombinedOutputLayerPrediction extends Prediction {
    /**
     * Genotype called by the model.
     */
    public String predictedGenotype;
    /**
     * Genotype called by the model.
     */
    public String trueGenotype;
    /**
     * The probability of the called genotype according to the model. Forecast probability.
     */
    public double overallProbability;

    public void inspectRecord(BaseInformationRecords.BaseInformation record) {
        trueGenotype=record.getTrueGenotype();
    }
}
