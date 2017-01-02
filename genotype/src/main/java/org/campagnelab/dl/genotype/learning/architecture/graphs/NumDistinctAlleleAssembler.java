package org.campagnelab.dl.genotype.learning.architecture.graphs;

import org.campagnelab.dl.framework.architecture.graphs.ComputationGraphAssembler;
import org.campagnelab.dl.framework.domains.DomainDescriptor;
import org.campagnelab.dl.framework.models.ModelPropertiesHelper;
import org.campagnelab.dl.framework.tools.TrainingArguments;
import org.campagnelab.dl.genotype.learning.GenotypeTrainingArguments;
import org.deeplearning4j.nn.api.OptimizationAlgorithm;
import org.deeplearning4j.nn.conf.ComputationGraphConfiguration;
import org.deeplearning4j.nn.conf.LearningRatePolicy;
import org.deeplearning4j.nn.conf.NeuralNetConfiguration;
import org.deeplearning4j.nn.conf.Updater;
import org.deeplearning4j.nn.conf.layers.DenseLayer;
import org.deeplearning4j.nn.conf.layers.OutputLayer;
import org.deeplearning4j.nn.graph.ComputationGraph;
import org.deeplearning4j.nn.weights.WeightInit;
import org.nd4j.linalg.lossfunctions.ILossFunction;

/**
 * Created by fac2003 on 12/20/16.
 */
public class NumDistinctAlleleAssembler implements ComputationGraphAssembler {
    private int numHiddenNodes;
    private LearningRatePolicy learningRatePolicy;
    private TrainingArguments arguments;
    private int numInputs;
    private String[] outputNames = new String[]{"numDistinctAlleles", "A", "T", "C", "G", "N", "I1", "I2", "I3", "I4", "I5", "metaData"};
    private int numLayers;

    private TrainingArguments args() {
        return arguments;
    }

   FeedForwardDenseLayerAssembler layerAssembler;

    @Override
    public void setArguments(TrainingArguments arguments) {
        this.arguments = arguments;
        this.numLayers = ((GenotypeTrainingArguments) arguments).numLayers;
        layerAssembler=new FeedForwardDenseLayerAssembler(arguments);
    }

    @Override
    public ComputationGraph createComputationalGraph(DomainDescriptor domainDescriptor) {
        layerAssembler.setLearningRatePolicy(learningRatePolicy);
        int numInputs = domainDescriptor.getNumInputs("input")[0];
        int numHiddenNodes = domainDescriptor.getNumHiddenNodes("firstDense");
        ComputationGraphConfiguration.GraphBuilder build = layerAssembler.assemble(numInputs, numHiddenNodes, numLayers);
        int numIn = layerAssembler.getNumOutputs();
        WeightInit WEIGHT_INIT = WeightInit.XAVIER;
        String lastDenseLayerName = layerAssembler.lastLayerName();

        build.addLayer("numDistinctAlleles", new OutputLayer.Builder(
                domainDescriptor.getOutputLoss("numDistinctAlleles"))
                .weightInit(WEIGHT_INIT)
                .activation("softmax").weightInit(WEIGHT_INIT).learningRateDecayPolicy(learningRatePolicy)
                .nIn(numIn)
                .nOut(domainDescriptor.getNumOutputs("numDistinctAlleles")[0]).build(), lastDenseLayerName).addInputs();
        build.addLayer("metaData", new OutputLayer.Builder(
                domainDescriptor.getOutputLoss("metaData"))
                .weightInit(WEIGHT_INIT)
                .activation("softmax").weightInit(WEIGHT_INIT).learningRateDecayPolicy(learningRatePolicy)
                .nIn(numIn)
                .nOut(
                        domainDescriptor.getNumOutputs("metaData")[0]
                ).build(), lastDenseLayerName);
        for (int i = 1; i < outputNames.length; i++) {
            build.addLayer(outputNames[i], new OutputLayer.Builder(
                    domainDescriptor.getOutputLoss(outputNames[i]))
                    .weightInit(WEIGHT_INIT)
                    .activation("softmax").weightInit(WEIGHT_INIT).learningRateDecayPolicy(learningRatePolicy)
                    .nIn(numIn).nOut(domainDescriptor.getNumOutputs(outputNames[i])[0]).build(), lastDenseLayerName);
        }
        ComputationGraphConfiguration conf = build
                .setOutputs(outputNames)
                .pretrain(false).backprop(true).build();

        return new ComputationGraph(conf);
    }

    @Override
    public void setNumInputs(String inputName, int... dimension) {

    }

    @Override
    public void setNumOutputs(String outputName, int... dimension) {

    }

    @Override
    public void setNumHiddenNodes(String componentName, int numHiddenNodes) {

    }

    @Override
    public String[] getInputNames() {
        return new String[]{"input"};

    }

    @Override
    public String[] getOutputNames() {
        return outputNames;
    }

    @Override
    public String[] getComponentNames() {
        return new String[]{"firstDense"};
    }

    @Override
    public void setLossFunction(String outputName, ILossFunction lossFunction) {
    }

    @Override
    public void saveProperties(ModelPropertiesHelper helper) {
    }
}
