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
    private final String[] hiddenLayerNames;
    private static final String[] lstmInputNames = new String[]{"from", "G1", "G2", "G3"};
    private final String[] inputVertexNames;
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
        hiddenLayerNames = new String[args().numLayers + args().numLSTMLayers];
        for (int i = 1; i <= args().numLayers; i++) {
            hiddenLayerNames[i] = "dense" + i;
        }
        int j = 0;
        for (int i = args().numLayers + 1; i < hiddenLayerNames.length; i++, j++) {
            hiddenLayerNames[i] = "lstmHidden" + j;
        }
        inputVertexNames = new String[lstmInputNames.length];
        for (int i = 0; i < lstmInputNames.length; i++) {
            inputVertexNames[i] = "lstm" + lstmInputNames[i].toUpperCase() + "OutputVertex";
        }
    }

    private GenotypeTrainingArguments args() {
        return arguments;
    }

    @Override
    public void setArguments(TrainingArguments arguments) {
        this.arguments = ((GenotypeTrainingArguments) arguments);
    }

    @Override
    public ComputationGraph createComputationalGraph(DomainDescriptor domainDescriptor) {
        LearningRatePolicy learningRatePolicy = LearningRatePolicy.Poly;
        int numInputs = domainDescriptor.getNumInputs("input")[0];
        int numLSTMInputs = domainDescriptor.getNumInputs("from")[0];
        int numHiddenNodes = domainDescriptor.getNumHiddenNodes("firstDense");
        WeightInit WEIGHT_INIT = WeightInit.XAVIER;
        if (numHiddenNodes >= 0) {
            throw new RuntimeException("model capacity is too small. At least some hidden nodes must be created.");
        }
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
                .weightInit(WeightInit.XAVIER).graphBuilder()
                .addInputs("input", "from", "G1", "G2", "G3")
                .setInputTypes(InputType.feedForward(numInputs),
                        InputType.recurrent(numLSTMInputs),
                        InputType.recurrent(numLSTMInputs),
                        InputType.recurrent(numLSTMInputs),
                        InputType.recurrent(numLSTMInputs));
        build.addLayer("lstmInput", new GravesLSTM.Builder()
                .nIn(numLSTMInputs)
                .nOut(numHiddenNodes)
                .activation("softsign")
                .build(), lstmInputNames);
        for (int i = 0; i < args().numLSTMLayers; i++) {
            String lstmPrevious = i == 0 ? "lstmInput" : "lstmHidden" + (i - 1);
            build.addLayer("lstmHidden" + i, new GravesLSTM.Builder()
                    .nIn(numHiddenNodes)
                    .nOut(numHiddenNodes)
                    .activation("softsign")
                    .build(), lstmPrevious);
        }
        build.addLayer("lstmOutput", new RnnOutputLayer.Builder(new LossMCXENT())
                .nIn(numHiddenNodes)
                .nOut(numInputs)
                .activation("softmax")
                .build(), "lstmHidden" + (args().numLSTMLayers - 1));
        for (int i = 0; i < lstmInputNames.length; i++) {
            build.addVertex(inputVertexNames[i], new LastTimeStepVertex(lstmInputNames[i]), "lstmOutput");
        }
        epsilon = 0.1;
        int numIn = numInputs;
        for (int i = 1; i <= args().numLayers; i++) {
            String lastDenseLayerName = "dense" + i;
            String previousLayerName = i == 1 ? "input" : "dense" + (i - 1);
            build.addLayer(lastDenseLayerName, new DenseLayer.Builder().nIn(numIn).nOut(numHiddenNodes)
                    .weightInit(WEIGHT_INIT)
                    .activation("relu").learningRateDecayPolicy(learningRatePolicy).epsilon(epsilon)
                    .build(), previousLayerName);
            numIn = numHiddenNodes;
        }
        String lastDenseLayerName = "dense" + args().numLayers;

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
