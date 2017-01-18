package org.campagnelab.dl.genotype.learning.architecture.graphs;

import org.campagnelab.dl.framework.tools.TrainingArguments;
import org.deeplearning4j.nn.api.OptimizationAlgorithm;
import org.deeplearning4j.nn.conf.ComputationGraphConfiguration;
import org.deeplearning4j.nn.conf.LearningRatePolicy;
import org.deeplearning4j.nn.conf.NeuralNetConfiguration;
import org.deeplearning4j.nn.conf.Updater;
import org.deeplearning4j.nn.conf.inputs.InputType;
import org.deeplearning4j.nn.conf.layers.DenseLayer;
import org.deeplearning4j.nn.weights.WeightInit;

/**
 * A helper class to assemble a configurable number of dense feed forward layers.
 * Created by fac2003 on 1/2/17.
 */
public class FeedForwardDenseLayerAssembler {
    private final ComputationGraphConfiguration.GraphBuilder build;
    private final int numGraphInputs;
    private TrainingArguments args;
    private LearningRatePolicy learningRatePolicy = LearningRatePolicy.Poly;
    private int numOutputs;
    private String lastLayerName;
    private static final double BUILDER_EPSILON = 1e-08d;
    private static final double LAYER_EPSILON = 0.1;
    private static final WeightInit WEIGHT_INIT = WeightInit.XAVIER;
    private static final float REDUCTION = 1f;


    public FeedForwardDenseLayerAssembler(TrainingArguments args, String... inputNames) {
        this.numGraphInputs = inputNames.length;
        setArgs(args);
        NeuralNetConfiguration.Builder graphBuilder = new NeuralNetConfiguration.Builder()
                .seed(args().seed)
                .iterations(1)
                .optimizationAlgo(OptimizationAlgorithm.STOCHASTIC_GRADIENT_DESCENT)
                .learningRate(args().learningRate)
                .updater(Updater.ADAGRAD)
                .epsilon(BUILDER_EPSILON)
                .lrPolicyDecayRate(0.5)
                .weightInit(WeightInit.XAVIER);
        if (args().regularizationRate != null) {
            graphBuilder.l2(args().regularizationRate);
        }
        if (args().dropoutRate != null) {
            graphBuilder.dropOut(args().dropoutRate);
            graphBuilder.setUseDropConnect(true);
        }
        build = graphBuilder.graphBuilder().addInputs(inputNames);
    }

    public FeedForwardDenseLayerAssembler(TrainingArguments args, ComputationGraphConfiguration.GraphBuilder build) {
        setArgs(args);
        this.build = build;
        this.numGraphInputs = build.getNetworkInputs().size();
    }

    private void setArgs(TrainingArguments args) {
        this.args = args;
    }

    private TrainingArguments args() {
        return args;
    }

    public ComputationGraphConfiguration.GraphBuilder assemble(int numInputs, int numHiddenNodes, int numLayers) {
        return assemble(numInputs, numHiddenNodes, numLayers, "input", 1);
    }

    public ComputationGraphConfiguration.GraphBuilder assemble(int numInputs, int numHiddenNodes, int numLayers,
                                                               String baseLayer, int startingIndex) {
        assert numHiddenNodes > 0 : "model capacity is too small. At least some hidden nodes must be created.";
        int minimum = (int) (numHiddenNodes * Math.pow(REDUCTION, 4));
        assert minimum > 2 : "Too much reduction, not enough outputs: ";
        int numIn = numInputs;
        int numOut = numHiddenNodes;
        String previousLayerName;
        String lastDenseLayerName = "no layers";
        for (int i = startingIndex; i < numLayers + startingIndex; i++) {
            numOut = numHiddenNodes;
            //     System.out.printf("layer %d numIn=%d numOut=%d%n", i, numIn, numOut);
            lastDenseLayerName = "dense" + i;
            previousLayerName = i == startingIndex ? baseLayer : "dense" + (i - 1);
            build.addLayer(lastDenseLayerName, new DenseLayer.Builder().nIn(numIn).nOut(numOut)
                    .weightInit(WEIGHT_INIT)
                    .activation("relu").learningRateDecayPolicy(learningRatePolicy).epsilon(LAYER_EPSILON)
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

    public void setInputTypes(InputType... inputTypes) {
        build.setInputTypes(inputTypes);
    }

    public ComputationGraphConfiguration.GraphBuilder getBuild() {
        return build;
    }

    public String lastLayerName() {
        return lastLayerName;
    }
}
