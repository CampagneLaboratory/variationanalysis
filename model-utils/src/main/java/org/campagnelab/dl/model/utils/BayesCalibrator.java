package org.campagnelab.dl.model.utils;

import java.io.IOException;

/**
 * Created by rct66 on 7/19/16.
 */
public class BayesCalibrator extends CalcCalibrator {

    public final double PRIOR_MUT_RATE = (double) 1/(4*(2e6));

    public BayesCalibrator(String modelPath, String prefix, boolean loadStats) throws IOException, ClassNotFoundException {
        super(modelPath, prefix, loadStats);
    }

    @Override
    public float calibrateProb(float prob) {
        double mutGreater = plantedMutSet.subSet(prob,1.1f).size() + 1;
        double unMutGreater = unMutSet.subSet(prob,1.1f).size() + 1;
        double  pMGreater = mutGreater/plantedMutSet.size();
        double  pUGreater = unMutGreater/ unMutSet.size();
        double bayes = pMGreater*PRIOR_MUT_RATE/pUGreater;
        return (float) bayes;
    }


}
