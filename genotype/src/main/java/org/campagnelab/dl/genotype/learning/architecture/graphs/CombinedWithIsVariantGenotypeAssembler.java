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
 * Created by fac2003 on 12/23/16.
 */
public class CombinedWithIsVariantGenotypeAssembler implements ComputationGraphAssembler {

    private int numHiddenNodes;
    private LearningRatePolicy learningRatePolicy;
    private TrainingArguments arguments;
    private int numInputs;
    private int numLayers;

    private TrainingArguments args() {
        return arguments;
    }


    @Override
    public void setArguments(TrainingArguments arguments) {
        this.arguments = arguments;
        this.numLayers = ((GenotypeTrainingArguments) arguments).numLayers;
    }

    @Override
    public ComputationGraph createComputationalGraph(DomainDescriptor domainDescriptor) {
        int numInputs = domainDescriptor.getNumInputs("input")[0];
        int numHiddenNodes = domainDescriptor.getNumHiddenNodes("firstDense");

        WeightInit WEIGHT_INIT = WeightInit.XAVIER;
        learningRatePolicy = LearningRatePolicy.Poly;
        float reduction = 1f;
        int minimum = (int) (numHiddenNodes * Math.pow(reduction, 4));
        assert minimum > 2 : "Too much reduction, not enough outputs: ";
        ComputationGraphConfiguration confBuilder = null;
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
                .weightInit(WeightInit.XAVIER).graphBuilder().addInputs("input");
        int numIn = numInputs;
        int numOut = numHiddenNodes;

        String lastDenseLayerName = "no layers";
        String previousLayerName = "input";
        epsilon = 0.1;
        for (int i = 1; i <= numLayers; i++) {
            numOut = numHiddenNodes;
            System.out.printf("layer %d numIn=%d numOut=%d%n", i, numIn, numOut);
            lastDenseLayerName = "dense" + i;

            previousLayerName = i == 1 ? "input" : "dense" + (i - 1);
            build.addLayer(lastDenseLayerName, new DenseLayer.Builder().nIn(numIn).nOut(numOut)
                    .weightInit(WEIGHT_INIT)
                    .activation("relu").learningRateDecayPolicy(learningRatePolicy).epsilon(epsilon)
                    .build(), previousLayerName);
            numIn = numOut;

        }

        build.addLayer("combined", new OutputLayer.Builder(
                domainDescriptor.getOutputLoss("combined"))
                .weightInit(WEIGHT_INIT)
                .activation("softmax").weightInit(WEIGHT_INIT).learningRateDecayPolicy(learningRatePolicy)
                .nIn(numIn)
                .nOut(domainDescriptor.getNumOutputs("combined")[0]).epsilon(epsilon).build(), lastDenseLayerName);
        build.addLayer("metaData", new OutputLayer.Builder(
                domainDescriptor.getOutputLoss("metaData"))
                .weightInit(WEIGHT_INIT).epsilon(epsilon)
                .activation("softmax").weightInit(WEIGHT_INIT).learningRateDecayPolicy(learningRatePolicy)
                .nIn(numIn)
                .nOut(
                        domainDescriptor.getNumOutputs("metaData")[0]
                ).epsilon(epsilon).build(), lastDenseLayerName);
        build.addLayer("isVariant", new OutputLayer.Builder(
                domainDescriptor.getOutputLoss("isVariant"))
                .weightInit(WEIGHT_INIT).epsilon(epsilon)
                .activation("softmax").weightInit(WEIGHT_INIT).learningRateDecayPolicy(learningRatePolicy)
                .nIn(numIn)
                .nOut(
                        domainDescriptor.getNumOutputs("isVariant")[0]
                ).epsilon(epsilon).build(), lastDenseLayerName);
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
        return new String[]{"combined", "metaData", "isVariant"};
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
