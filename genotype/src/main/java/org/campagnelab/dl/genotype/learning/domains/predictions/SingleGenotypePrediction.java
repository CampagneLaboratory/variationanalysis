package org.campagnelab.dl.genotype.learning.domains.predictions;

import org.campagnelab.dl.framework.domains.prediction.Prediction;
import org.campagnelab.dl.varanalysis.protobuf.BaseInformationRecords;

/**
 * Created by rct66 on 11/12/16.
 */
public class SingleGenotypePrediction extends Prediction {


    public String predictedSingleGenotype;
    public double probabilityIsCalled;
    public boolean trueIsCalled;


}

