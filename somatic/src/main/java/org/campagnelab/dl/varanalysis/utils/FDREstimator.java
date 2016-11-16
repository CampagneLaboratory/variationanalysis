package org.campagnelab.dl.varanalysis.utils;

import java.io.IOException;

/**
 * Created by rct66 on 7/19/16.
 */
public class FDREstimator extends CalcCalibrator {

    public FDREstimator(String modelPath, String prefix, boolean loadStats) throws IOException, ClassNotFoundException {
        super(modelPath, prefix, loadStats);
    }

    @Override
    public float calibrateProb(float prob) {
        double unMutGreater = unMutSet.subSet(prob,1.1f).size();
        double unMutTotal = unMutSet.size();
        double FDR = unMutGreater/(unMutTotal);
        return (float) FDR;
    }


}
