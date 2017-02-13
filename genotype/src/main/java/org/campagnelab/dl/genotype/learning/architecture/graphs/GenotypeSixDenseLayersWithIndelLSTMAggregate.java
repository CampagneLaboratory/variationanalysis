package org.campagnelab.dl.genotype.learning.architecture.graphs;

import org.campagnelab.dl.framework.architecture.graphs.ComputationGraphAssembler;
import org.campagnelab.dl.framework.domains.DomainDescriptor;
import org.campagnelab.dl.framework.models.ModelPropertiesHelper;
import org.campagnelab.dl.framework.tools.TrainingArguments;
import org.campagnelab.dl.genotype.learning.GenotypeTrainingArguments;
import org.deeplearning4j.nn.conf.ComputationGraphConfiguration;
import org.deeplearning4j.nn.conf.LearningRatePolicy;
import org.deeplearning4j.nn.conf.graph.MergeVertex;
import org.deeplearning4j.nn.conf.graph.rnn.LastTimeStepVertex;
import org.deeplearning4j.nn.conf.inputs.InputType;
import org.deeplearning4j.nn.conf.layers.GravesLSTM;
import org.deeplearning4j.nn.conf.layers.OutputLayer;
import org.deeplearning4j.nn.graph.ComputationGraph;
import org.deeplearning4j.nn.weights.WeightInit;
import org.nd4j.linalg.lossfunctions.ILossFunction;

/**
 * Created by joshuacohen on 2/6/17.
 */
public class GenotypeSixDenseLayersWithIndelLSTMAggregate extends GenotypeAssembler implements ComputationGraphAssembler {
    private final String[] outputNames;
    private final String[] inputNames;
    private static final String[] lstmInputNames = new String[]{"indel"};
    private static final WeightInit WEIGHT_INIT = WeightInit.XAVIER;
    private static final LearningRatePolicy LEARNING_RATE_POLICY = LearningRatePolicy.Poly;
    private static final GenotypeSixDenseLayersWithIndelLSTMAggregate.OutputType DEFAULT_OUTPUT_TYPE = GenotypeSixDenseLayersWithIndelLSTMAggregate.OutputType.DISTINCT_ALLELES;
    private final String combined;
    private final GenotypeSixDenseLayersWithIndelLSTMAggregate.OutputType outputType;
    private final boolean addTrueGenotypeLabels;

    public enum OutputType {
        HOMOZYGOUS,
        DISTINCT_ALLELES,
        COMBINED,
    }

    private GenotypeTrainingArguments arguments;

    public GenotypeSixDenseLayersWithIndelLSTMAggregate() {
        this(DEFAULT_OUTPUT_TYPE, false, false, false);
    }

    public GenotypeSixDenseLayersWithIndelLSTMAggregate(GenotypeSixDenseLayersWithIndelLSTMAggregate.OutputType outputType, boolean hasIsVariant) {
        this(outputType, hasIsVariant, false, false);
    }

    public GenotypeSixDenseLayersWithIndelLSTMAggregate(GenotypeSixDenseLayersWithIndelLSTMAggregate.OutputType outputType,
                                                        boolean hasIsVariant, boolean fixRef, boolean addTrueGenotypeLabels) {
        this.outputType = outputType;
        this.hasIsVariant = hasIsVariant;
        this.addTrueGenotypeLabels = addTrueGenotypeLabels;
        combined = fixRef ? "combinedRef" : "combined";
        switch (outputType) {
            case DISTINCT_ALLELES:
                if (hasIsVariant) {
                    if (addTrueGenotypeLabels) {
                        outputNames = new String[]{"numDistinctAlleles", "A", "T", "C", "G", "N",
                                "I1", "I2", "I3", "I4", "I5", "metaData", "isVariant", "trueGenotype"};
                    } else {
                        outputNames = new String[]{"numDistinctAlleles", "A", "T", "C", "G", "N",
                                "I1", "I2", "I3", "I4", "I5", "metaData", "isVariant"};
                    }
                } else {
                    if (addTrueGenotypeLabels) {
                        outputNames = new String[]{"numDistinctAlleles", "A", "T", "C", "G", "N",
                                "I1", "I2", "I3", "I4", "I5", "metaData", "trueGenotype"};
                    } else {
                        outputNames = new String[]{"numDistinctAlleles", "A", "T", "C", "G", "N",
                                "I1", "I2", "I3", "I4", "I5", "metaData"};
                    }
                }
                break;
            case COMBINED:
                if (hasIsVariant) {
                    if (addTrueGenotypeLabels) {
                        outputNames = new String[]{combined, "metaData", "isVariant", "trueGenotype"};
                    } else {
                        outputNames = new String[]{combined, "metaData", "isVariant"};
                    }
                } else {
                    if (addTrueGenotypeLabels) {
                        outputNames = new String[]{combined, "metaData", "trueGenotype"};
                    } else {
                        outputNames = new String[]{combined, "metaData"};
                    }
                }
                break;
            case HOMOZYGOUS:
                if (hasIsVariant) {
                    if (addTrueGenotypeLabels) {
                        outputNames = new String[]{"homozygous", "A", "T", "C", "G", "N",
                                "I1", "I2", "I3", "I4", "I5", "metaData", "isVariant", "trueGenotype"};
                    } else {
                        outputNames = new String[]{"homozygous", "A", "T", "C", "G", "N",
                                "I1", "I2", "I3", "I4", "I5", "metaData", "isVariant"};
                    }
                } else {
                    if (addTrueGenotypeLabels) {
                        outputNames = new String[]{"homozygous", "A", "T", "C", "G", "N",
                                "I1", "I2", "I3", "I4", "I5", "metaData", "trueGenotype"};
                    } else {
                        outputNames = new String[]{"homozygous", "A", "T", "C", "G", "N",
                                "I1", "I2", "I3", "I4", "I5", "metaData"};
                    }
                }
                break;
            default:
                throw new RuntimeException("Output type not recognized");
        }
        if (addTrueGenotypeLabels) {
            inputNames = new String[]{"input", "indel", "trueGenotypeInput"};
        } else {
            inputNames = new String[]{"input", "indel"};
        }
    }

    private GenotypeTrainingArguments args() {
        return arguments;
    }

    @Override
    public void setArguments(TrainingArguments arguments) {
        this.arguments = ((GenotypeTrainingArguments) arguments);
    }

    private void addOutputLayers(ComputationGraphConfiguration.GraphBuilder build, DomainDescriptor domainDescriptor,
                                 String lastDenseLayerName, int numIn, int numLSTMLayers, int numLSTMHiddenNodes,
                                 int numLSTMInputs) {
        if (outputType == GenotypeSixDenseLayersWithIndelLSTMAggregate.OutputType.HOMOZYGOUS || outputType == GenotypeSixDenseLayersWithIndelLSTMAggregate.OutputType.DISTINCT_ALLELES) {
            if (outputType == GenotypeSixDenseLayersWithIndelLSTMAggregate.OutputType.DISTINCT_ALLELES) {
                build.addLayer("numDistinctAlleles", new OutputLayer.Builder(domainDescriptor.getOutputLoss("homozygous"))
                        .weightInit(WEIGHT_INIT)
                        .activation("softmax").weightInit(WEIGHT_INIT)
                        .learningRateDecayPolicy(LEARNING_RATE_POLICY)
                        .nIn(numIn)
                        .nOut(domainDescriptor.getNumOutputs("numDistinctAlleles")[0])
                        .build(), lastDenseLayerName);
            } else {
                build.addLayer("homozygous", new OutputLayer.Builder(domainDescriptor.getOutputLoss("homozygous"))
                        .weightInit(WEIGHT_INIT)
                        .activation("softmax").weightInit(WEIGHT_INIT).learningRateDecayPolicy(LEARNING_RATE_POLICY)
                        .nIn(numIn).nOut(11).build(), lastDenseLayerName);
            }
            int endingIndex = hasIsVariant ? outputNames.length - 2 : outputNames.length - 1;
            for (int i = 1; i <= endingIndex; i++) {
                build.addLayer(outputNames[i], new OutputLayer.Builder(domainDescriptor.getOutputLoss(outputNames[i]))
                        .weightInit(WEIGHT_INIT)
                        .activation("softmax").weightInit(WEIGHT_INIT).learningRateDecayPolicy(LEARNING_RATE_POLICY)
                        .nIn(numIn).nOut(2).build(), lastDenseLayerName);
            }
        } else if (outputType == GenotypeSixDenseLayersWithIndelLSTMAggregate.OutputType.COMBINED) {
            build.addLayer(combined, new OutputLayer.Builder(domainDescriptor.getOutputLoss(combined))
                    .weightInit(WEIGHT_INIT)
                    .activation(combined).weightInit(WEIGHT_INIT).learningRateDecayPolicy(LEARNING_RATE_POLICY)
                    .nIn(numIn)
                    .nOut(domainDescriptor.getNumOutputs(combined)[0]).build(), lastDenseLayerName);
        }
        appendMetaDataLayer(domainDescriptor, LEARNING_RATE_POLICY, build, numIn, WEIGHT_INIT, lastDenseLayerName);
        appendIsVariantLayer(domainDescriptor, LEARNING_RATE_POLICY, build, numIn, WEIGHT_INIT, lastDenseLayerName);
        appendTrueGenotypeLayers(build, addTrueGenotypeLabels, lastDenseLayerName, domainDescriptor, WEIGHT_INIT, LEARNING_RATE_POLICY,
                numLSTMLayers, numIn, numLSTMHiddenNodes);
    }

    @Override
    public ComputationGraph createComputationalGraph(DomainDescriptor domainDescriptor) {
        int numInputs = domainDescriptor.getNumInputs("input")[0];
        int numLSTMInputs = domainDescriptor.getNumInputs("indel")[0];
        int numHiddenNodes = domainDescriptor.getNumHiddenNodes("firstDense");
        int numLSTMHiddenNodes = domainDescriptor.getNumHiddenNodes("lstmLayer");
        int numLSTMLayers = Math.max(1, args().numLSTMLayers);
        int numLayers = Math.max(1, args().numLayers);
        FeedForwardDenseLayerAssembler assembler = new FeedForwardDenseLayerAssembler(args());
        assembler.setLearningRatePolicy(LEARNING_RATE_POLICY);
        assembler.initializeBuilder(getInputNames());
        assembler.setInputTypes(getInputTypes(domainDescriptor));
        ComputationGraphConfiguration.GraphBuilder build = assembler.getBuild();
        String lstmLayerName = "no layer";
        for (int i = 0; i < numLSTMLayers; i++) {
            lstmLayerName = "lstmindel_" + i;
            String lstmPreviousLayerName = i == 0 ? "indel" : "lstmindel_" + (i - 1);
            int numLSTMInputNodes = i == 0 ? numLSTMInputs : numLSTMHiddenNodes;
            build.addLayer(lstmLayerName, new GravesLSTM.Builder()
                    .nIn(numLSTMInputNodes)
                    .nOut(numLSTMHiddenNodes)
                    .build(), lstmPreviousLayerName);
        }
        build.addVertex("lstmindelLastTimeStepVertex", new LastTimeStepVertex("indel"),
                lstmLayerName);
        String[] mergeInputs = new String[lstmInputNames.length + 1];
        for (int i = 0; i < lstmInputNames.length; i++) {
            mergeInputs[i] = "lstm" + lstmInputNames[i] + "LastTimeStepVertex";
        }
        mergeInputs[lstmInputNames.length] = "input";
        build.addVertex("lstmFeedForwardMerge", new MergeVertex(), mergeInputs);
        int numInputsToDenseAfterMerge = numInputs + (lstmInputNames.length  * numLSTMHiddenNodes);
        assembler.assemble(numInputsToDenseAfterMerge, numHiddenNodes,
                numLayers, "lstmFeedForwardMerge", 1);
        String lastDenseLayerName = assembler.lastLayerName();
        int numIn = assembler.getNumOutputs();
        addOutputLayers(build, domainDescriptor, lastDenseLayerName, numIn, numLSTMLayers, numLSTMHiddenNodes,
                numLSTMInputs);
        ComputationGraphConfiguration conf = build
                .setOutputs(outputNames)
                .build();
        System.out.println(conf);
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
        return inputNames;
    }

    private InputType[] getInputTypes(DomainDescriptor domainDescriptor) {
        String[] inputNames = getInputNames();
        InputType[] inputTypes = new InputType[inputNames.length];
        for (int i = 0; i < inputNames.length; i++) {
            switch(inputNames[i]) {
                case "input":
                    inputTypes[i] = InputType.feedForward(domainDescriptor.getNumInputs(inputNames[i])[0]);
                    break;
                case "indel":
                    inputTypes[i] = InputType.recurrent(domainDescriptor.getNumInputs(inputNames[i])[0]);
                    break;
                default:
                    throw new RuntimeException("Invalid input to computation graph");
            }
        }
        return inputTypes;
    }

    @Override
    public String[] getOutputNames() {
        return outputNames;
    }

    @Override
    public String[] getComponentNames() {
        return new String[]{"firstDense", "lstmLayer"};
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
