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
 * Created by fac2003 on 12/23/16.
 */
public class CombinedGenotypeAssembler extends GenotypeAssembler implements ComputationGraphAssembler {


    private TrainingArguments arguments;
    private int numLayers;
    private FeedForwardDenseLayerAssembler layerAssembler;
    private boolean fixRef = false;

    public CombinedGenotypeAssembler() {
        this(false);
    }

    public CombinedGenotypeAssembler(boolean hasIsVariant) {
        this.hasIsVariant = hasIsVariant;
    }

    public CombinedGenotypeAssembler(boolean hasIsVariant, boolean fixRef) {
        this.hasIsVariant = hasIsVariant;
        this.fixRef = fixRef;
    }

    private TrainingArguments args() {
        return arguments;
    }


    @Override
    public void setArguments(TrainingArguments arguments) {
        this.arguments = arguments;
        this.numLayers = ((GenotypeTrainingArguments) arguments).numLayers;
        layerAssembler = new FeedForwardDenseLayerAssembler(arguments, getInputNames());
    }

    @Override
    public ComputationGraph createComputationalGraph(DomainDescriptor domainDescriptor) {
        LearningRatePolicy learningRatePolicy = LearningRatePolicy.Poly;
        layerAssembler.setLearningRatePolicy(learningRatePolicy);
        int numInputs = domainDescriptor.getNumInputs("input")[0];
        int numHiddenNodes = domainDescriptor.getNumHiddenNodes("firstDense");

        ComputationGraphConfiguration.GraphBuilder build = layerAssembler.assemble(numInputs,
                numHiddenNodes, numLayers);

        int numIn = layerAssembler.getNumOutputs();
        WeightInit WEIGHT_INIT = WeightInit.XAVIER;
        String lastDenseLayerName = layerAssembler.lastLayerName();
        String combined = fixRef?"combinedRef":"combined";
        build.addLayer(combined, new OutputLayer.Builder(
                domainDescriptor.getOutputLoss(combined))
                .weightInit(WEIGHT_INIT)
                .activation(combined).weightInit(WEIGHT_INIT).learningRateDecayPolicy(learningRatePolicy)
                .nIn(numIn)
                .nOut(domainDescriptor.getNumOutputs(combined)[0]).build(), lastDenseLayerName);
        appendMetaDataLayer(domainDescriptor, learningRatePolicy, build, numIn, WEIGHT_INIT, lastDenseLayerName);
        appendIsVariantLayer(domainDescriptor, learningRatePolicy, build, numIn, WEIGHT_INIT, lastDenseLayerName);
        ComputationGraphConfiguration conf = build
                .setOutputs(getOutputNames())
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
        String combined =  fixRef ? "combinedRef" : "combined";

        if (hasIsVariant) {
            return new String[]{combined, "metaData", "isVariant"};
        } else {
            return new String[]{combined, "metaData"};
        }
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
        helper.put(this.getClass().getCanonicalName() + ".numLayers", numLayers);
    }
}
