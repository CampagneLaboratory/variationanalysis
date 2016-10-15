package org.campagnelab.dl.varanalysis.learning.architecture;

import org.deeplearning4j.nn.api.OptimizationAlgorithm;
import org.deeplearning4j.nn.conf.LearningRatePolicy;
import org.deeplearning4j.nn.conf.MultiLayerConfiguration;
import org.deeplearning4j.nn.conf.NeuralNetConfiguration;
import org.deeplearning4j.nn.conf.Updater;
import org.deeplearning4j.nn.conf.layers.DenseLayer;
import org.deeplearning4j.nn.conf.layers.OutputLayer;

/**
 * A network with six dense layers. This is the second neural net architecture we tried for detecting
 * somatic variations.
 * <p>
 * Created by fac2003 on 6/10/16.
 *
 * @author Fabien Campagne
 */
public class SixDenseLayersForGPU extends AbstractNeuralNetAssembler implements NeuralNetAssembler {


    public MultiLayerConfiguration createNetwork() {
        learningRatePolicy = LearningRatePolicy.Score;

        NeuralNetConfiguration.ListBuilder confBuilder = null;
        NeuralNetConfiguration.Builder netBuilder = new NeuralNetConfiguration.Builder()
                .seed(seed)
                .iterations(1)
                .optimizationAlgo(OptimizationAlgorithm.STOCHASTIC_GRADIENT_DESCENT)
                .learningRate(learningRate).regularization(regularization).l2(regularizationRate)
                .updater(Updater.ADAGRAD);

        if (dropOut) {
            netBuilder.dropOut(dropOutRate);
            netBuilder.setUseDropConnect(true);
        }
        numHiddenNodes = 0;
        // at least as many hidden nodes as input nodes, but multiple number of 64.
        while (numInputs > numHiddenNodes) {
            numHiddenNodes += 64;
        }
        confBuilder = netBuilder.list()
                .layer(0, new DenseLayer.Builder().nIn(numInputs).nOut(numHiddenNodes)
                        .weightInit(WEIGHT_INIT)
                        .activation("relu").learningRateDecayPolicy(learningRatePolicy)
                        .build())
                .layer(1, new DenseLayer.Builder().nIn(numHiddenNodes).nOut((int) (numHiddenNodes))
                        .weightInit(WEIGHT_INIT)
                        .activation("relu").learningRateDecayPolicy(learningRatePolicy)
                        .build())
                .layer(2, new DenseLayer.Builder().nIn((int) (numHiddenNodes)).nOut((int) (numHiddenNodes))
                        .weightInit(WEIGHT_INIT).learningRateDecayPolicy(learningRatePolicy)
                        .activation("relu")
                        .build())
                .layer(3, new DenseLayer.Builder().nIn((int) (numHiddenNodes)).nOut((int) (numHiddenNodes))
                        .weightInit(WEIGHT_INIT).learningRateDecayPolicy(learningRatePolicy)
                        .activation("relu")
                        .build())
                .layer(4, new DenseLayer.Builder().nIn((int) (numHiddenNodes)).nOut((int) (numHiddenNodes))
                        .weightInit(WEIGHT_INIT).learningRateDecayPolicy(learningRatePolicy)
                        .activation("relu")
                        .build())
                .layer(5, new OutputLayer.Builder(lossFunction)
                        .weightInit(WEIGHT_INIT)
                        .activation("softmax").weightInit(WEIGHT_INIT).learningRateDecayPolicy(learningRatePolicy)
                        .nIn((int) (numHiddenNodes)).nOut(numOutputs).build())
                .pretrain(false).backprop(true);
        return confBuilder.build();

    }
}
