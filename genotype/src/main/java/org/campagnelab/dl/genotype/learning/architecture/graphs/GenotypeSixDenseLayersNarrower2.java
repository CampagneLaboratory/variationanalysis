package org.campagnelab.dl.genotype.learning.architecture.graphs;

import org.campagnelab.dl.framework.architecture.graphs.ComputationGraphAssembler;
import org.campagnelab.dl.framework.domains.DomainDescriptor;
import org.campagnelab.dl.framework.models.ModelPropertiesHelper;
import org.campagnelab.dl.framework.tools.TrainingArguments;
import org.campagnelab.dl.genotype.learning.GenotypeTrainingArguments;
import org.deeplearning4j.nn.conf.ComputationGraphConfiguration;
import org.deeplearning4j.nn.conf.LearningRatePolicy;
import org.deeplearning4j.nn.conf.WorkspaceMode;
import org.deeplearning4j.nn.conf.layers.OutputLayer;
import org.deeplearning4j.nn.graph.ComputationGraph;
import org.deeplearning4j.nn.weights.WeightInit;
import org.nd4j.linalg.lossfunctions.ILossFunction;

/**
 * A computational graph with six dense layers and two outputs: probabilityIsCalled of somatic mutation and frequency of the
 * mutation (0 when no mutation).
 * <p>
 *
 * @author Remi Torracinta
 */
public class GenotypeSixDenseLayersNarrower2 extends GenotypeAssembler implements ComputationGraphAssembler {


    private int numHiddenNodes;
    private TrainingArguments arguments;
    private int numInputs;
    private String[] outputNames;
    private FeedForwardDenseLayerAssembler layerAssembler;
    private int numLayers;

    private TrainingArguments args() {
        return arguments;
    }

    public GenotypeSixDenseLayersNarrower2() {
        this(false);
    }

    public GenotypeSixDenseLayersNarrower2(boolean hasIsVariant) {
        this.hasIsVariant = hasIsVariant;
        if (hasIsVariant) {
            outputNames = new String[]{"homozygous", "A", "T", "C", "G", "N",
                    "I1", "I2", "I3", "I4", "I5", "metaData"};
        } else {
            outputNames = new String[]{"homozygous", "A", "T", "C", "G", "N",
                    "I1", "I2", "I3", "I4", "I5", "metaData", "isVariant"};
        }
    }

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

        build.addLayer("homozygous", new OutputLayer.Builder(
                domainDescriptor.getOutputLoss("homozygous"))
                .weightInit(WEIGHT_INIT)
                .activation("softmax").weightInit(WEIGHT_INIT).learningRateDecayPolicy(learningRatePolicy)
                .nIn(numIn).nOut(11).build(), lastDenseLayerName);
        for (int i = 1; i < outputNames.length; i++) {
            build.addLayer(outputNames[i], new OutputLayer.Builder(
                    domainDescriptor.getOutputLoss(outputNames[i]))
                    .weightInit(WEIGHT_INIT)
                    .activation("softmax").weightInit(WEIGHT_INIT).learningRateDecayPolicy(learningRatePolicy)
                    .nIn(numIn).nOut(2).build(), lastDenseLayerName);
        }
        appendMetaDataLayer(domainDescriptor, learningRatePolicy, build, numIn, WEIGHT_INIT, lastDenseLayerName);
        appendIsVariantLayer(domainDescriptor, learningRatePolicy, build, numIn, WEIGHT_INIT, lastDenseLayerName);

        ComputationGraphConfiguration conf = build
                .setOutputs(outputNames)
                .build();
        conf.setTrainingWorkspaceMode(WorkspaceMode.SINGLE);
        conf.setInferenceWorkspaceMode(WorkspaceMode.SINGLE);
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
