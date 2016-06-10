package org.campagnelab.dl.varanalysis.learning.architecture;

import org.deeplearning4j.nn.api.OptimizationAlgorithm;
import org.deeplearning4j.nn.conf.MultiLayerConfiguration;
import org.deeplearning4j.nn.conf.NeuralNetConfiguration;
import org.deeplearning4j.nn.conf.Updater;
import org.deeplearning4j.nn.conf.layers.DenseLayer;
import org.deeplearning4j.nn.conf.layers.OutputLayer;
import org.deeplearning4j.nn.weights.WeightInit;
import org.nd4j.linalg.lossfunctions.LossFunctions;

/**
 * A network with six dense layers. This is the second neural net architecture we tried for detecting
 * somatic variations.
 * <p>
 * Created by fac2003 on 6/10/16.
 *
 * @author Fabien Campagne
 */
public class SixDenseLayers extends AbstractNeuralNetAssembler implements NeuralNetAssembler {



    public MultiLayerConfiguration createNetwork() {
        MultiLayerConfiguration conf = new NeuralNetConfiguration.Builder()
                .seed(seed)
                .iterations(10)
                .optimizationAlgo(OptimizationAlgorithm.STOCHASTIC_GRADIENT_DESCENT)
                .learningRate(learningRate).regularization(regularization).l2(regularizationRate)
                .updater(Updater.ADAGRAD)
                .list()
                .layer(0, new DenseLayer.Builder().nIn(numInputs).nOut(numHiddenNodes)
                        .weightInit(WeightInit.XAVIER)
                        .activation("relu")
                        .build())
                .layer(1, new DenseLayer.Builder().nIn(numHiddenNodes).nOut(numHiddenNodes)
                        .weightInit(WeightInit.XAVIER)
                        .activation("relu")
                        .build())
                .layer(2, new DenseLayer.Builder().nIn(numHiddenNodes).nOut(numHiddenNodes)
                        .weightInit(WeightInit.XAVIER)
                        .activation("relu")
                        .build())
                .layer(3, new DenseLayer.Builder().nIn(numHiddenNodes).nOut(numHiddenNodes)
                        .weightInit(WeightInit.XAVIER)
                        .activation("relu")
                        .build())
                .layer(4, new DenseLayer.Builder().nIn(numHiddenNodes).nOut(numHiddenNodes)
                        .weightInit(WeightInit.XAVIER)
                        .activation("relu")
                        .build())
                .layer(5, new OutputLayer.Builder(LossFunctions.LossFunction.NEGATIVELOGLIKELIHOOD)
                        .weightInit(WeightInit.XAVIER)
                        .activation("softmax").weightInit(WeightInit.XAVIER)
                        .nIn(numHiddenNodes).nOut(numOutputs).build())
                .pretrain(false).backprop(true).build();
        return conf;

    }
}
