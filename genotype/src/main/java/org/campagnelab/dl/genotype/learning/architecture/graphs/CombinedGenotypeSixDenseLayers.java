package org.campagnelab.dl.genotype.learning.architecture.graphs;

import org.campagnelab.dl.framework.architecture.graphs.ComputationGraphAssembler;
import org.campagnelab.dl.framework.domains.DomainDescriptor;
import org.campagnelab.dl.framework.models.ModelPropertiesHelper;
import org.campagnelab.dl.framework.tools.TrainingArguments;
import org.deeplearning4j.nn.api.OptimizationAlgorithm;
import org.deeplearning4j.nn.conf.ComputationGraphConfiguration;
import org.deeplearning4j.nn.conf.LearningRatePolicy;
import org.deeplearning4j.nn.conf.NeuralNetConfiguration;
import org.deeplearning4j.nn.conf.Updater;
import org.deeplearning4j.nn.conf.layers.DenseLayer;
import org.deeplearning4j.nn.conf.layers.OutputLayer;
import org.deeplearning4j.nn.graph.ComputationGraph;
import org.deeplearning4j.nn.weights.WeightInit;
import org.nd4j.linalg.lossfunctions.LossFunctions;

/**
 * A computational graph with six dense layers and two outputs: probabilityIsCalled of somatic mutation and frequency of the
 * mutation (0 when no mutation).
 * <p>
 *
 * @author Remi Torracinta
 */
public class CombinedGenotypeSixDenseLayers implements ComputationGraphAssembler {


    private int numHiddenNodes;
    private LearningRatePolicy learningRatePolicy;
    private TrainingArguments arguments;
    private int numInputs;
    private String[] outputNames = new String[]{"combined"};
    private TrainingArguments args() {
        return arguments;
    };

    @Override
    public void setArguments(TrainingArguments arguments) {
        this.arguments = arguments;
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
        NeuralNetConfiguration.Builder graphBuilder = new NeuralNetConfiguration.Builder()
                .seed(args().seed)
                .iterations(1)
                .optimizationAlgo(OptimizationAlgorithm.STOCHASTIC_GRADIENT_DESCENT)
                .learningRate(args().learningRate)
                .updater(Updater.ADAGRAD);
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
                .weightInit(WeightInit.XAVIER).graphBuilder().addInputs("input")
                .addLayer("dense1", new DenseLayer.Builder().nIn(numInputs).nOut(numHiddenNodes)
                        .weightInit(WEIGHT_INIT)
                        .activation("relu").learningRateDecayPolicy(learningRatePolicy)
                        .build(), "input")
                .addLayer("dense2", new DenseLayer.Builder().nIn(numHiddenNodes).nOut((int) (numHiddenNodes * reduction))
                        .weightInit(WEIGHT_INIT)
                        .activation("relu").learningRateDecayPolicy(learningRatePolicy)
                        .build(), "dense1")
                .addLayer("dense3", new DenseLayer.Builder().nIn((int) (numHiddenNodes * reduction)).nOut((int) (numHiddenNodes * Math.pow(reduction, 2)))
                        .weightInit(WEIGHT_INIT).learningRateDecayPolicy(learningRatePolicy)
                        .activation("relu")
                        .build(), "dense2")
                .addLayer("dense4", new DenseLayer.Builder().nIn((int) (numHiddenNodes * Math.pow(reduction, 2))).nOut((int) (numHiddenNodes * Math.pow(reduction, 3)))
                        .weightInit(WEIGHT_INIT).learningRateDecayPolicy(learningRatePolicy)
                        .activation("relu")
                        .build(), "dense3")
                .addLayer("dense5", new DenseLayer.Builder().nIn((int) (numHiddenNodes * Math.pow(reduction, 3))).nOut((int) (numHiddenNodes * Math.pow(reduction, 4)))
                        .weightInit(WEIGHT_INIT).learningRateDecayPolicy(learningRatePolicy)
                        .activation("relu")
                        .build(), "dense4");
        build.addLayer("combined", new OutputLayer.Builder(
                domainDescriptor.getOutputLoss("combined"))
                .weightInit(WEIGHT_INIT)
                .activation("softmax").weightInit(WEIGHT_INIT).learningRateDecayPolicy(learningRatePolicy)
                .nIn((int) (numHiddenNodes * Math.pow(reduction, 4))).nOut(4).build(), "dense5");
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
    public void setLossFunction(String outputName, LossFunctions.LossFunction lossFunction) {
    }

    @Override
    public void saveProperties(ModelPropertiesHelper helper) {
    }
}
