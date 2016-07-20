package org.campagnelab.dl.model.utils;

import org.campagnelab.dl.model.utils.CalcCalibrator;

import java.io.IOException;
import java.net.MalformedURLException;

/**
 * Created by rct66 on 7/19/16.
 */
public class FDREstimator extends CalcCalibrator {

    public FDREstimator(String modelPath, boolean loadStats) throws IOException, ClassNotFoundException {
        super(modelPath, loadStats);
    }

    @Override
    public float calibrateProb(float prob) {
        double mutGreater = plantedMutSet.subSet(prob,1.1f).size() + 1;
        double unMutGreater = unMutSet.subSet(prob,1.1f).size() + 1;

        double FDR = unMutGreater/(unMutGreater+mutGreater);
        return (float) FDR;
    }


}
