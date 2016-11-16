package org.campagnelab.dl.somatic.learning.architecture;

import org.campagnelab.dl.framework.architecture.nets.NeuralNetAssembler;
import org.deeplearning4j.nn.api.OptimizationAlgorithm;
import org.deeplearning4j.nn.conf.LearningRatePolicy;
import org.deeplearning4j.nn.conf.MultiLayerConfiguration;
import org.deeplearning4j.nn.conf.NeuralNetConfiguration;
import org.deeplearning4j.nn.conf.Updater;
import org.deeplearning4j.nn.conf.layers.BatchNormalization;
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
public class SixDenseLayersNarrower2BatchNormalization extends AbstractNeuralNetAssembler implements NeuralNetAssembler {




    public MultiLayerConfiguration createNetwork() {
        learningRatePolicy = LearningRatePolicy.Poly;
        float reduction = 0.65f;
        int minimum = (int) (numHiddenNodes * Math.pow(reduction, 4));
        assert minimum > numOutputs : "Too much reduction, not enough outputs: ";
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

        confBuilder = netBuilder.list()
                .layer(0, new DenseLayer.Builder().nIn(numInputs).nOut(numHiddenNodes)
                        .weightInit(WEIGHT_INIT)
                        .activation("relu").learningRateDecayPolicy(learningRatePolicy)
                        .build())
                .layer(1,new BatchNormalization.Builder().nOut( (numHiddenNodes))
                        .build())
                .layer(2, new DenseLayer.Builder().nIn(numHiddenNodes).nOut((int) (numHiddenNodes * reduction))
                        .weightInit(WEIGHT_INIT)
                        .activation("relu").learningRateDecayPolicy(learningRatePolicy)
                        .build())
                .layer(3,new BatchNormalization.Builder().nOut((int) (numHiddenNodes * reduction))
                        .build())
                .layer(4, new DenseLayer.Builder().nIn((int) (numHiddenNodes * reduction)).nOut((int) (numHiddenNodes * Math.pow(reduction, 2)))
                        .weightInit(WEIGHT_INIT).learningRateDecayPolicy(learningRatePolicy)
                        .activation("relu")
                        .build())
                .layer(5, new DenseLayer.Builder().nIn((int) (numHiddenNodes * Math.pow(reduction, 2))).nOut((int) (numHiddenNodes * Math.pow(reduction, 3)))
                        .weightInit(WEIGHT_INIT).learningRateDecayPolicy(learningRatePolicy)
                        .activation("relu")
                        .build())

                .layer(6, new DenseLayer.Builder().nIn((int) (numHiddenNodes * Math.pow(reduction, 3))).nOut((int) (numHiddenNodes * Math.pow(reduction, 4)))
                        .weightInit(WEIGHT_INIT).learningRateDecayPolicy(learningRatePolicy)
                        .activation("relu")
                        .build())

                .layer(7, new OutputLayer.Builder(lossFunction)
                        .weightInit(WEIGHT_INIT)
                        .activation("softmax").weightInit(WEIGHT_INIT).learningRateDecayPolicy(learningRatePolicy)
                        .nIn((int) (numHiddenNodes * Math.pow(reduction, 4))).nOut(numOutputs).build())
                .pretrain(false).backprop(true);
        return confBuilder.build();

    }
}
