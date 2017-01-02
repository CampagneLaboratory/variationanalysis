package org.campagnelab.dl.genotype.learning.architecture.graphs;

import org.campagnelab.dl.framework.tools.TrainingArguments;
import org.deeplearning4j.nn.api.OptimizationAlgorithm;
import org.deeplearning4j.nn.conf.ComputationGraphConfiguration;
import org.deeplearning4j.nn.conf.LearningRatePolicy;
import org.deeplearning4j.nn.conf.NeuralNetConfiguration;
import org.deeplearning4j.nn.conf.Updater;
import org.deeplearning4j.nn.conf.layers.DenseLayer;
import org.deeplearning4j.nn.weights.WeightInit;

/**
 * A helper class to assemble a configurable number of dense feed forward layers.
 * Created by fac2003 on 1/2/17.
 */
public class FeedForwardDenseLayerAssembler {
    private TrainingArguments args;
    private LearningRatePolicy learningRatePolicy;
    private int numOutputs;
    private String lastLayerName;

    public FeedForwardDenseLayerAssembler(TrainingArguments args) {
        this.args = args;
    }

    TrainingArguments args() {
        return args;
    }

    ComputationGraphConfiguration.GraphBuilder assemble(int numInputs, int numHiddenNodes, int numLayers) {
        assert numHiddenNodes > 0 : "model capacity is too small. At least some hidden nodes must be created.";
        WeightInit WEIGHT_INIT = WeightInit.XAVIER;
        float reduction = 1f;
        int minimum = (int) (numHiddenNodes * Math.pow(reduction, 4));
        assert minimum > 2 : "Too much reduction, not enough outputs: ";
        ComputationGraphConfiguration confBuilder = null;
        double epsilon = 1e-08d;
        NeuralNetConfiguration.Builder graphBuilder = new NeuralNetConfiguration.Builder()
                .seed(args().seed)
                .iterations(1)
                .optimizationAlgo(OptimizationAlgorithm.STOCHASTIC_GRADIENT_DESCENT)
                .learningRate(args().learningRate)
                .updater(Updater.ADAGRAD);
        graphBuilder.epsilon(epsilon);
        if (args().regularizationRate != null) {
            graphBuilder.l2(args().regularizationRate);
        }
        if (args().dropoutRate != null) {
            graphBuilder.dropOut(args().dropoutRate);
            graphBuilder.setUseDropConnect(true);
        }
        NeuralNetConfiguration.Builder graphConfiguration = graphBuilder.lrPolicyDecayRate(0.5)
                .optimizationAlgo(OptimizationAlgorithm.STOCHASTIC_GRADIENT_DESCENT).iterations(1)
                .learningRate(args().learningRate)
                .seed(args().seed);
        if (args().regularizationRate != null) {
            graphConfiguration.regularization(args().regularizationRate != null);
        }
        if (args().dropoutRate != null) {
            graphConfiguration.dropOut(args().dropoutRate);
            graphConfiguration.setUseDropConnect(true);
        }
        ComputationGraphConfiguration.GraphBuilder build = graphConfiguration
                .weightInit(WeightInit.XAVIER).graphBuilder().addInputs("input");
        int numIn = numInputs;
        int numOut = numHiddenNodes;

        String lastDenseLayerName = "no layers";
        String previousLayerName = "input";
        epsilon = 0.1;
        for (int i = 1; i <= numLayers; i++) {
            numOut = numHiddenNodes;
            //     System.out.printf("layer %d numIn=%d numOut=%d%n", i, numIn, numOut);
            lastDenseLayerName = "dense" + i;

            previousLayerName = i == 1 ? "input" : "dense" + (i - 1);
            build.addLayer(lastDenseLayerName, new DenseLayer.Builder().nIn(numIn).nOut(numOut)
                    .weightInit(WEIGHT_INIT)
                    .activation("relu").learningRateDecayPolicy(learningRatePolicy).epsilon(epsilon)
                    .build(), previousLayerName);
            numIn = numOut;

        }

        this.numOutputs = numOut;
        this.lastLayerName = lastDenseLayerName;
        return build;

    }

    public int getNumOutputs() {
        return numOutputs;
    }

    public void setLearningRatePolicy(LearningRatePolicy learningRatePolicy) {
        this.learningRatePolicy = learningRatePolicy;
    }


    public String lastLayerName() {
        return lastLayerName;
    }
}
