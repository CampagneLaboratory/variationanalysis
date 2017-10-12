package org.campagnelab.dl.genotype.learning.architecture.graphs;

import it.unimi.dsi.fastutil.objects.ObjectArrayList;
import org.campagnelab.dl.framework.architecture.graphs.ComputationGraphAssembler;
import org.campagnelab.dl.framework.domains.DomainDescriptor;
import org.campagnelab.dl.framework.mappers.MappedDimensions;
import org.campagnelab.dl.framework.models.ModelPropertiesHelper;
import org.campagnelab.dl.framework.tools.TrainingArguments;
import org.campagnelab.dl.genotype.learning.GenotypeTrainingArguments;
import org.campagnelab.dl.genotype.performance.BEDHelper;
import org.deeplearning4j.nn.conf.ComputationGraphConfiguration;
import org.deeplearning4j.nn.conf.LearningRatePolicy;
import org.deeplearning4j.nn.conf.graph.MergeVertex;
import org.deeplearning4j.nn.conf.graph.rnn.LastTimeStepVertex;
import org.deeplearning4j.nn.conf.inputs.InputType;
import org.deeplearning4j.nn.conf.layers.GravesLSTM;
import org.deeplearning4j.nn.conf.layers.LSTM;
import org.deeplearning4j.nn.conf.layers.OutputLayer;
import org.deeplearning4j.nn.graph.ComputationGraph;
import org.deeplearning4j.nn.weights.WeightInit;
import org.nd4j.linalg.activations.impl.ActivationSoftmax;
import org.nd4j.linalg.lossfunctions.ILossFunction;

import java.util.Set;

/**
 * A computational graph with an LSTM over bases in a genomic segment [input is 2D].
 * <p>
 *
 * @author Fabien Campagne
 */
public class GenotypeSegmentsLSTM extends GenotypeAssembler implements ComputationGraphAssembler {

    private final String[] outputNames;
    private static final WeightInit WEIGHT_INIT = WeightInit.XAVIER;
    private static final LearningRatePolicy LEARNING_RATE_POLICY = LearningRatePolicy.Poly;
    private static final GenotypeSixDenseLayersWithIndelLSTM.OutputType DEFAULT_OUTPUT_TYPE = GenotypeSixDenseLayersWithIndelLSTM.OutputType.DISTINCT_ALLELES;
    private TrainingArguments arguments;
    private int numLayers;
    private FeedForwardDenseLayerAssembler layerAssembler;
    LearningRatePolicy learningRatePolicy = LearningRatePolicy.Poly;
    private ObjectArrayList<String> componentNames = new ObjectArrayList<>();

    public GenotypeSegmentsLSTM() {

        outputNames = new String[]{"genotype"};

    }

    private GenotypeTrainingArguments args() {
        return (GenotypeTrainingArguments) arguments;
    }


    @Override
    public void setArguments(TrainingArguments arguments) {
        this.arguments = arguments;
        this.numLayers = ((GenotypeTrainingArguments) arguments).numLayers;
        layerAssembler = new FeedForwardDenseLayerAssembler(arguments);

    }


    @Override
    public ComputationGraph createComputationalGraph(DomainDescriptor domainDescriptor) {
        // the number of elements per time-step (# mapped features for each base):
        final int[] numInputs = domainDescriptor.getNumInputs("input");
        int numLSTMInputs = numInputs[1];
        System.out.printf("GenotypeSegmentsLSTM getNumInputs: sequence-length=%d, num-float-per-base=%d%n",
                numInputs[0], numInputs[1]);
        int numHiddenNodes = domainDescriptor.getNumHiddenNodes("lstm");
        int numLSTMLayers = Math.max(1, args().numLSTMLayers);

        FeedForwardDenseLayerAssembler assembler = new FeedForwardDenseLayerAssembler(args());
        assembler.setLearningRatePolicy(LEARNING_RATE_POLICY);
        assembler.initializeBuilder(getInputNames());
        assembler.setInputTypes(getInputTypes(domainDescriptor));
        ComputationGraphConfiguration.GraphBuilder build = assembler.getBuild();
        String lstmInputName = "input";
        String lstmLayerName = "no layer";
        for (int i = 0; i < numLSTMLayers; i++) {
            lstmLayerName = "lstm_" + lstmInputName + "_" + i;
            String lstmPreviousLayerName = i == 0 ? lstmInputName : "lstm" + lstmInputName + "_" + (i - 1);
            int numLSTMInputNodes = i == 0 ? numLSTMInputs : numHiddenNodes;
            build.addLayer(lstmLayerName, new LSTM.Builder()
                    .nIn(numLSTMInputNodes)
                    .nOut(numHiddenNodes)
                    .weightInit(WEIGHT_INIT)
                    .build(), lstmPreviousLayerName);
            componentNames.add(lstmLayerName);
        }

        String lastDenseLayerName = assembler.lastLayerName();
        int numIn = assembler.getNumOutputs();
        build.addLayer("softmaxGenotype", new OutputLayer.Builder(
                domainDescriptor.getOutputLoss("Genotype"))
                .weightInit(WEIGHT_INIT)
                .activation(new ActivationSoftmax()).weightInit(WEIGHT_INIT).learningRateDecayPolicy(learningRatePolicy)
                .nIn(numIn)
                .nOut(domainDescriptor.getNumOutputs("Genotype")[0]).build(), lastDenseLayerName).addInputs();

        ComputationGraphConfiguration conf = build
                .setOutputs(outputNames)
                .build();
        // System.out.println(conf);
        return new ComputationGraph(conf);
    }

    private InputType getInputTypes(DomainDescriptor domainDescriptor) {
        final MappedDimensions inputDimensions = domainDescriptor.getFeatureMapper("input").dimensions();
        System.out.printf("GenotypeSegmentsLSTM dimensions: sequence-length=%d, num-float-per-base=%d%n",
                inputDimensions.numElements(0), inputDimensions.numElements());
        return InputType.recurrent(inputDimensions.numElements(0), inputDimensions.numElements());
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
        return componentNames.toArray(new String[componentNames.size()]);
    }

    @Override
    public void setLossFunction(String outputName, ILossFunction lossFunction) {
    }

    @Override
    public void saveProperties(ModelPropertiesHelper helper) {
    }
}
