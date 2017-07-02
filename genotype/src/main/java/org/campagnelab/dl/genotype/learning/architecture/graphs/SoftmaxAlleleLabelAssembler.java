package org.campagnelab.dl.genotype.learning.architecture.graphs;

import org.campagnelab.dl.framework.architecture.graphs.ComputationGraphAssembler;
import org.campagnelab.dl.framework.domains.DomainDescriptor;
import org.campagnelab.dl.framework.models.ModelPropertiesHelper;
import org.campagnelab.dl.framework.tools.TrainingArguments;
import org.campagnelab.dl.genotype.learning.GenotypeTrainingArguments;
import org.deeplearning4j.nn.conf.ComputationGraphConfiguration;
import org.deeplearning4j.nn.conf.LearningRatePolicy;
import org.deeplearning4j.nn.conf.layers.OutputLayer;
import org.deeplearning4j.nn.graph.ComputationGraph;
import org.deeplearning4j.nn.weights.WeightInit;
import org.nd4j.linalg.lossfunctions.ILossFunction;

/**
 * A computation graph with combined softmax allele predictions. Predicts the combination of alleles that are called.
 * This assembler stores the labels in metaData.
 * Created by fac2003 on 2/7/16.
 */
public class SoftmaxAlleleLabelAssembler extends GenotypeAssembler implements ComputationGraphAssembler {

    private TrainingArguments arguments;
    private String[] outputNames;
    private int numLayers;

    public SoftmaxAlleleLabelAssembler() {
        this(false);
    }

    public SoftmaxAlleleLabelAssembler(boolean hasIsVariant) {
        this.hasIsVariant = hasIsVariant;
        if (hasIsVariant) {
            outputNames = new String[]{"softmaxGenotype",  "metaData", "isVariant"};
        } else {
            outputNames = new String[]{"softmaxGenotype",  "metaData",};
        }
    }

    private TrainingArguments args() {
        return arguments;
    }

    private FeedForwardDenseLayerAssembler layerAssembler;

    @Override
    public void setArguments(TrainingArguments arguments) {
        this.arguments = arguments;
        this.numLayers = ((GenotypeTrainingArguments) arguments).numLayers;
        layerAssembler = new FeedForwardDenseLayerAssembler(arguments);
    }

    @Override
    public ComputationGraph createComputationalGraph(DomainDescriptor domainDescriptor) {
        LearningRatePolicy learningRatePolicy = LearningRatePolicy.Poly;
        layerAssembler.setLearningRatePolicy(learningRatePolicy);
        layerAssembler.initializeBuilder();
        int numInputs = domainDescriptor.getNumInputs("input")[0];
        int numHiddenNodes = domainDescriptor.getNumHiddenNodes("firstDense");
        ComputationGraphConfiguration.GraphBuilder build = layerAssembler.assemble(numInputs, numHiddenNodes, numLayers);
        int numIn = layerAssembler.getNumOutputs();
        WeightInit WEIGHT_INIT = WeightInit.XAVIER;
        String lastDenseLayerName = layerAssembler.lastLayerName();

        build.addLayer("softmaxGenotype", new OutputLayer.Builder(
                domainDescriptor.getOutputLoss("softmaxGenotype"))
                .weightInit(WEIGHT_INIT)
                .activation("softmax").weightInit(WEIGHT_INIT).learningRateDecayPolicy(learningRatePolicy)
                .nIn(numIn)
                .nOut(domainDescriptor.getNumOutputs("softmaxGenotype")[0]).build(), lastDenseLayerName).addInputs();


        appendMetaDataLayer(domainDescriptor, learningRatePolicy, build, numIn, WEIGHT_INIT, lastDenseLayerName);
        appendIsVariantLayer(domainDescriptor, learningRatePolicy, build, numIn, WEIGHT_INIT, lastDenseLayerName);

        ComputationGraphConfiguration conf = build
                .setOutputs(outputNames)
                .build();

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
