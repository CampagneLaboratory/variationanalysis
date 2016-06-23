package org.campagnelab.dl.varanalysis.learning.architecture;

import org.deeplearning4j.earlystopping.EarlyStoppingConfiguration;
import org.deeplearning4j.earlystopping.EarlyStoppingResult;
import org.deeplearning4j.earlystopping.listener.EarlyStoppingListener;
import org.deeplearning4j.nn.multilayer.MultiLayerNetwork;

/**
 * Created by rct66 on 6/21/16.
 */
public class EStatusListener implements EarlyStoppingListener<org.deeplearning4j.nn.multilayer.MultiLayerNetwork> {
    @Override
    public void onStart(EarlyStoppingConfiguration<MultiLayerNetwork> esConfig, MultiLayerNetwork net) {
        System.out.println("Training has started on:" + net.toString());
    }

    @Override
    public void onEpoch(int epochNum, double score, EarlyStoppingConfiguration<MultiLayerNetwork> esConfig, MultiLayerNetwork net) {
        System.out.println("Epoch" + Integer.toString(epochNum) + "has completed with score:" + Double.toString(score));
    }

    @Override
    public void onCompletion(EarlyStoppingResult<MultiLayerNetwork> esResult) {
        System.out.println("Training has completed with result:" + esResult.toString());
    }
}
