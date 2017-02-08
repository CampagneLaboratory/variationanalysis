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
 * A computation graph with individual allele predictions (one for each base and indel), and the number
 * of distinct allele to call. An optional isVariant output can be included (see constructor flag).
 * This assembler stores the labels in metaData.
 * Created by fac2003 on 12/20/16.
 */
public class NumDistinctAlleleAssembler extends GenotypeAssembler implements ComputationGraphAssembler {

    private TrainingArguments arguments;
    private String[] outputNames;
    private int numLayers;

    public NumDistinctAlleleAssembler() {
        this(false);
    }

    public NumDistinctAlleleAssembler(boolean hasIsVariant) {
        this.hasIsVariant = hasIsVariant;
        if (hasIsVariant) {
            outputNames = new String[]{"numDistinctAlleles", "A", "T", "C", "G", "N",
                    "I1", "I2", "I3", "I4", "I5", "metaData", "isVariant"};
        } else {
            outputNames = new String[]{"numDistinctAlleles", "A", "T", "C", "G", "N",
                    "I1", "I2", "I3", "I4", "I5", "metaData"};
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

        build.addLayer("numDistinctAlleles", new OutputLayer.Builder(
                domainDescriptor.getOutputLoss("numDistinctAlleles"))
                .weightInit(WEIGHT_INIT)
                .activation("softmax").weightInit(WEIGHT_INIT).learningRateDecayPolicy(learningRatePolicy)
                .nIn(numIn)
                .nOut(domainDescriptor.getNumOutputs("numDistinctAlleles")[0]).build(), lastDenseLayerName).addInputs();

        for (int i = 1; i < outputNames.length; i++) {
            build.addLayer(outputNames[i], new OutputLayer.Builder(
                    domainDescriptor.getOutputLoss(outputNames[i]))
                    .weightInit(WEIGHT_INIT)
                    .activation("softmax").weightInit(WEIGHT_INIT).learningRateDecayPolicy(learningRatePolicy)
                    .nIn(numIn).nOut(domainDescriptor.getNumOutputs(outputNames[i])[0]).build(), lastDenseLayerName);
        }
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
