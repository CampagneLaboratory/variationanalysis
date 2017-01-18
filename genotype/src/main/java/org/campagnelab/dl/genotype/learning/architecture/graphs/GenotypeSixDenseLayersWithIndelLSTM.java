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
import org.deeplearning4j.nn.conf.graph.MergeVertex;
import org.deeplearning4j.nn.conf.graph.rnn.LastTimeStepVertex;
import org.deeplearning4j.nn.conf.inputs.InputType;
import org.deeplearning4j.nn.conf.layers.DenseLayer;
import org.deeplearning4j.nn.conf.layers.GravesLSTM;
import org.deeplearning4j.nn.conf.layers.OutputLayer;
import org.deeplearning4j.nn.conf.layers.RnnOutputLayer;
import org.deeplearning4j.nn.graph.ComputationGraph;
import org.deeplearning4j.nn.weights.WeightInit;
import org.nd4j.linalg.lossfunctions.ILossFunction;
import org.nd4j.linalg.lossfunctions.impl.LossMCXENT;

/**
 * Created by joshuacohen on 1/12/17.
 */
public class GenotypeSixDenseLayersWithIndelLSTM extends GenotypeAssembler implements ComputationGraphAssembler {
    private final String[] outputNames;
    private String[] hiddenLayerNames;
    private static final String[] lstmInputNames = new String[]{"from", "G1", "G2", "G3"};
    private static final WeightInit WEIGHT_INIT = WeightInit.XAVIER;
    private static final LearningRatePolicy LEARNING_RATE_POLICY = LearningRatePolicy.Poly;

    private GenotypeTrainingArguments arguments;


    public GenotypeSixDenseLayersWithIndelLSTM() {
        this(false);
    }

    public GenotypeSixDenseLayersWithIndelLSTM(boolean hasIsVariant) {
        this.hasIsVariant = hasIsVariant;
        if (hasIsVariant) {
            outputNames = new String[]{"homozygous", "A", "T", "C", "G", "N",
                    "I1", "I2", "I3", "I4", "I5", "metaData"};
        } else {
            outputNames = new String[]{"homozygous", "A", "T", "C", "G", "N",
                    "I1", "I2", "I3", "I4", "I5", "metaData", "isVariant"};
        }
    }

    private GenotypeTrainingArguments args() {
        return arguments;
    }

    public void setHiddenLayerNames() {
        hiddenLayerNames = new String[args().numLayers + args().numPreVertexLayers + (args().numLSTMLayers * lstmInputNames.length)];
        for (int i = 0; i < args().numLayers + args().numPreVertexLayers; i++) {
            hiddenLayerNames[i] = "dense" + i;
        }
        int c = args().numLayers + args().numPreVertexLayers;
        for (int i = 0; i < args().numLSTMLayers; i++) {
            for (String lstmInputName : lstmInputNames) {
                hiddenLayerNames[c++] = "lstm" + lstmInputName + "Hidden" + i;
            }
        }
    }
    @Override
    public void setArguments(TrainingArguments arguments) {
        this.arguments = ((GenotypeTrainingArguments) arguments);
    }

    @Override
    public ComputationGraph createComputationalGraph(DomainDescriptor domainDescriptor) {
        int numInputs = domainDescriptor.getNumInputs("input")[0];
        int numLSTMInputs = domainDescriptor.getNumInputs("from")[0];
        int numHiddenNodes = domainDescriptor.getNumHiddenNodes("firstDense");
        if (hiddenLayerNames == null) {
            setHiddenLayerNames();
        }
        FeedForwardDenseLayerAssembler assembler = new FeedForwardDenseLayerAssembler(args(), "input", "from", "G1", "G2", "G3");
        assembler.setInputTypes(InputType.feedForward(numInputs), InputType.recurrent(numLSTMInputs), InputType.recurrent(numLSTMInputs),
                InputType.recurrent(numLSTMInputs), InputType.recurrent(numLSTMInputs));
        ComputationGraphConfiguration.GraphBuilder build = assembler.getBuild();
        for (String lstmInputName : lstmInputNames) {
            String lstmInputLayerName = "lstm" + lstmInputName + "Input";
            build.addLayer(lstmInputLayerName, new GravesLSTM.Builder()
                    .nIn(numLSTMInputs)
                    .nOut(numHiddenNodes)
                    .activation("softsign")
                    .build(), lstmInputName);
            for (int i = 0; i < args().numLSTMLayers; i++) {
                String lstmPrevious = i == 0 ? lstmInputLayerName :  "lstm" + lstmInputName + "Hidden" + (i - 1);
                build.addLayer("lstm" + lstmInputName + "Hidden" + i, new GravesLSTM.Builder()
                        .nIn(numHiddenNodes)
                        .nOut(numHiddenNodes)
                        .activation("softsign")
                        .build(), lstmPrevious);
            }
            String lstmOutputLayerName = "lstm" + lstmInputName + "Output";
            build.addLayer(lstmOutputLayerName, new RnnOutputLayer.Builder(new LossMCXENT())
                    .nIn(numHiddenNodes)
                    .nOut(args().numLSTMOutputs)
                    .activation("softmax")
                    .build(), "lstm" + lstmInputName + "Hidden" + (args().numLSTMLayers - 1));
            build.addVertex("lstm" + lstmInputName + "LastTimeStepVertex", new LastTimeStepVertex(lstmInputName),
                    lstmOutputLayerName);
        }
        assembler.assemble(numInputs, numHiddenNodes, args().numPreVertexLayers);
        String[] mergeInputs = new String[lstmInputNames.length + 1];
        for (int i = 0; i < lstmInputNames.length; i++) {
            mergeInputs[i] = "lstm" + lstmInputNames[i] + "LastTimeStepVertex";
        }
        mergeInputs[lstmInputNames.length] = assembler.lastLayerName();
        build.addVertex("lstmFeedForwardMerge", new MergeVertex(), mergeInputs);
        assembler.assemble(numHiddenNodes + args().numLSTMOutputs, numHiddenNodes,
                args().numLayers, "lstmFeedForwardMerge", args().numPreVertexLayers + 1);
        // TODO: check proper outputs
        String lastDenseLayerName = assembler.lastLayerName();
        int numIn = assembler.getNumOutputs();
        build.addLayer("homozygous", new OutputLayer.Builder(
                domainDescriptor.getOutputLoss("homozygous"))
                .weightInit(WEIGHT_INIT)
                .activation("softmax").weightInit(WEIGHT_INIT).learningRateDecayPolicy(LEARNING_RATE_POLICY)
                .nIn(numIn).nOut(11).build(), lastDenseLayerName);
        for (int i = 1; i < outputNames.length; i++) {
            build.addLayer(outputNames[i], new OutputLayer.Builder(
                    domainDescriptor.getOutputLoss(outputNames[i]))
                    .weightInit(WEIGHT_INIT)
                    .activation("softmax").weightInit(WEIGHT_INIT).learningRateDecayPolicy(LEARNING_RATE_POLICY)
                    .nIn(numIn).nOut(2).build(), lastDenseLayerName);
        }
        appendMetaDataLayer(domainDescriptor, LEARNING_RATE_POLICY, build, numIn, WEIGHT_INIT, lastDenseLayerName);
        appendIsVariantLayer(domainDescriptor, LEARNING_RATE_POLICY, build, numIn, WEIGHT_INIT, lastDenseLayerName);
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
        return new String[]{"input", "from", "G1", "G2", "G3"};
    }

    @Override
    public String[] getOutputNames() {
        return outputNames;
    }

    @Override
    public String[] getComponentNames() {
        return hiddenLayerNames;
    }

    @Override
    public void setLossFunction(String outputName, ILossFunction lossFunction) {

    }

    @Override
    public void saveProperties(ModelPropertiesHelper helper) {
        helper.put(this.getClass().getCanonicalName() + ".numLayers", args().numLayers);
        helper.put(this.getClass().getCanonicalName() + ".numLstmLayers", args().numLSTMLayers);
    }
}
