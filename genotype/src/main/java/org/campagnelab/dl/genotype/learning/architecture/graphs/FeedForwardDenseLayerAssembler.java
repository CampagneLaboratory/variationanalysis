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
    private TrainingArguments args;
    private LearningRatePolicy learningRatePolicy = LearningRatePolicy.Poly;
    private int numOutputs;
    private String lastLayerName;
    private ComputationGraphConfiguration.GraphBuilder build;
    private static final double BUILDER_EPSILON = 1e-08d;
    private static final double LAYER_EPSILON = 0.1;
    private static final WeightInit WEIGHT_INIT = WeightInit.XAVIER;
    private static final float REDUCTION = 1f;
    private double modelCapacity=1;
    private double reductionRate=1;

    public FeedForwardDenseLayerAssembler(TrainingArguments args) {
        this.args = args;
    }

    TrainingArguments args() {
        return args;
    }

    public void initializeBuilder(String... inputNames) {
        if (inputNames.length == 0) {
            inputNames = new String[]{"input"};
        }
        NeuralNetConfiguration.Builder graphBuilder = new NeuralNetConfiguration.Builder()
                .seed(args().seed)
                .iterations(1)
                .optimizationAlgo(OptimizationAlgorithm.STOCHASTIC_GRADIENT_DESCENT)
                .learningRate(args().learningRate)
                .updater(Updater.ADAGRAD)
                .epsilon(BUILDER_EPSILON)
                .lrPolicyDecayRate(0.5)
                .weightInit(WEIGHT_INIT);
        if (args().regularizationRate != null) {
            graphBuilder.l2(args().regularizationRate);
            graphBuilder.regularization(args().regularizationRate != null);
        }
        if (args().dropoutRate != null) {
            graphBuilder.dropOut(args().dropoutRate);
            graphBuilder.setUseDropConnect(true);
        }
        modelCapacity=args().modelCapacity;
        reductionRate=args().reductionRate;
        build = graphBuilder.graphBuilder().addInputs(inputNames);
    }

    public ComputationGraphConfiguration.GraphBuilder assemble(int numInputs, int numHiddenNodes, int numLayers) {
        return assemble(numInputs, numHiddenNodes, numLayers, "input", 1);
    }

    ComputationGraphConfiguration.GraphBuilder assemble(int numInputs, int numHiddenNodes, int numLayers, String baseLayer, int startingIndex) {
        assert numHiddenNodes > 0 : "model capacity is too small. At least some hidden nodes must be created.";
        WeightInit WEIGHT_INIT = WeightInit.XAVIER;
        float reduction = 1f;
        int minimum = (int) (numHiddenNodes * Math.pow(reduction, numLayers));
        assert minimum > 2 : "Too much reduction, not enough outputs: ";
        int numIn = numInputs;
        int numOut = (int) (numHiddenNodes * modelCapacity);
        String previousLayerName;
        String lastDenseLayerName = "no layers";
        for (int i = startingIndex; i < startingIndex + numLayers; i++) {

            //     System.out.printf("layer %d numIn=%d numOut=%d%n", i, numIn, numOut);
            lastDenseLayerName = "dense" + i;
            previousLayerName = i == startingIndex ? baseLayer : "dense" + (i - 1);
            numOut = (int) (numHiddenNodes * Math.pow(reductionRate, i) * modelCapacity);
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
