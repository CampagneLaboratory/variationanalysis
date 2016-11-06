package org.campagnelab.dl.varanalysis.learning;

import it.unimi.dsi.fastutil.floats.FloatArraySet;
import it.unimi.dsi.fastutil.floats.FloatSet;
import org.apache.commons.io.FileUtils;
import org.campagnelab.dl.model.utils.mappers.FeatureMapper;
import org.campagnelab.dl.model.utils.mappers.LabelMapper;
import org.campagnelab.dl.model.utils.models.ModelLoader;
import org.campagnelab.dl.varanalysis.learning.architecture.ComputationalGraphAssembler;
import org.campagnelab.dl.varanalysis.learning.models.ModelPropertiesHelper;
import org.campagnelab.dl.varanalysis.learning.models.PerformanceLogger;
import org.campagnelab.dl.varanalysis.tools.ConditionRecordingTool;
import org.campagnelab.goby.baseinfo.SequenceBaseInformationReader;
import org.deeplearning4j.earlystopping.EarlyStoppingResult;
import org.deeplearning4j.nn.api.Layer;
import org.deeplearning4j.nn.api.Model;
import org.deeplearning4j.nn.graph.ComputationGraph;
import org.nd4j.linalg.api.ndarray.INDArray;
import org.nd4j.linalg.dataset.api.iterator.MultiDataSetIterator;
import org.nd4j.linalg.lossfunctions.LossFunctions;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Date;
import java.util.Properties;
import java.util.function.Function;

/**
 * An abstract tool to train computational graphs. Implements early stopping. This class defines
 * several abstract methods that must be implemented to adapt training to different problems.
 */
public abstract class TrainModel extends ConditionRecordingTool<TrainingArguments> {

    static private Logger LOG = LoggerFactory.getLogger(TrainModel.class);

    private String directory;
    private double bestScore;
    private long time;

    protected abstract FeatureMapper getFeatureMapper(String inputName);

    protected abstract LabelMapper getLabelMapper(String outputName);

    protected abstract Function<String[], MultiDataSetIterator> getIteratorFunction();

    protected abstract ComputationalGraphAssembler getComputationalGraph();

    protected abstract int[] getNumInputs(String inputName);

    protected abstract int[] getNumOutputs(String outputName);

    protected abstract int getNumHiddenNodes(String componentName);

    protected abstract LossFunctions.LossFunction getOutputLoss(String outputName);

    protected abstract EarlyStoppingResult<ComputationGraph> train(ComputationGraph graph,
                                                                   MultiDataSetIterator async)
            throws IOException;

    protected PerformanceLogger performanceLogger;

    protected FeatureMapper featureMapper = null;
    protected LabelMapper labelMapper = null;

    @Override
    public void execute() {
        if (args().getTrainingSets().length == 0) {
            System.err.println("You must provide training datasets.");
        }

        try {
            featureMapper = getFeatureMapper("input");
            labelMapper = getLabelMapper("output");
            execute(featureMapper, args().getTrainingSets(), args().miniBatchSize);
        } catch (IOException e) {
            System.err.println("An exception occured. Details may be provided below");
            e.printStackTrace();
        }

    }

    public void execute(FeatureMapper featureCalculator, String trainingDataset[], int miniBatchSize) throws IOException {
        if (args().previousModelPath != null) {
            System.out.println(String.format("Resuming training with %s model parameters from %s %n", args().previousModelName, args().previousModelPath));
        }
        time = new Date().getTime();
        System.out.println("time: " + time);
        System.out.println("epochs: " + args().maxEpochs);
        System.out.println(featureCalculator.getClass().getTypeName());
        directory = "models/" + Long.toString(time);
        FileUtils.forceMkdir(new File(directory));

        // Assemble the training iterator:
        MultiDataSetIterator trainingIterator = getIteratorFunction().apply(trainingDataset);
        System.out.println("Estimating scaling parameters:");
        //Load the training data:


        ComputationalGraphAssembler assembler = getComputationalGraph();
        assembler.setArguments(args());
        for (String inputName : assembler.getInputNames()) {
            assembler.setNumInputs(inputName, getNumInputs(inputName));
        }
        for (String outputName : assembler.getOutputNames()) {
            assembler.setNumOutputs(outputName, getNumOutputs(outputName));
            assembler.setLossFunction(outputName, getOutputLoss(outputName));
        }
        for (String componentName : assembler.getOutputNames()) {
            assembler.setNumHiddenNodes(componentName, getNumHiddenNodes(componentName));
        }

        ComputationGraph net = assembler.createComputationalGraph();
        net.init();
        if (args().previousModelPath != null) {
            // Load the parameters of a previously trained model and set them on the new model to continue
            // training where we left it off. Note that models must have the same architecture or setting
            // parameters will fail.

            ModelLoader loader = new ModelLoader(args().previousModelPath);
            Model savedNetwork = loader.loadModel(args().previousModelName);
            ComputationGraph savedGraph = savedNetwork instanceof ComputationGraph ?
                    (ComputationGraph) savedNetwork :
                    null;
            if (savedNetwork == null || savedGraph.getUpdater() == null || savedGraph.params() == null) {
                System.err.println("Unable to load model or updater from " + args().previousModelPath);
            } else {
                net.setUpdater(savedGraph.getUpdater());
                net.setParams(savedNetwork.params());
            }
        }

        //Print the  number of parameters in the graph (and for each layer)
        Layer[] layers = net.getLayers();
        int totalNumParams = 0;
        for (int i = 0; i < layers.length; i++) {
            int nParams = layers[i].numParams();
            System.out.println("Number of parameters in layer " + i + ": " + nParams);
            totalNumParams += nParams;
        }
        System.out.println("Total number of network parameters: " + totalNumParams);

        writeProperties(assembler, this);
        performanceLogger = new PerformanceLogger(directory);
        EarlyStoppingResult<ComputationGraph> result = train(net, trainingIterator);

        //Print out the results:
        System.out.println("Termination reason: " + result.getTerminationReason());
        System.out.println("Termination details: " + result.getTerminationDetails());
        System.out.println("Total epochs: " + result.getTotalEpochs());
        System.out.println("Best epoch number: " + result.getBestModelEpoch());
        System.out.println("Score at best epoch: " + performanceLogger.getBestScore());
        System.out.println("AUC at best epoch: " + performanceLogger.getBestAUC());

        writeProperties(assembler, this);
        writeBestScoreFile();
        System.out.println("Model completed, saved at time: " + time);
        performanceLogger.write();
        resultValues().put("AUC", performanceLogger.getBestAUC());
        resultValues().put("score", performanceLogger.getBestScore());
        resultValues().put("bestModelEpoch", performanceLogger.getBestEpoch("bestAUC"));
        resultValues().put("model-time", time);
    }


    protected void writeBestScoreFile() throws IOException {

        FileWriter scoreWriter = new FileWriter(directory + "/bestScore");
        scoreWriter.append(Double.toString(bestScore));
        scoreWriter.close();
    }

    protected void writeProperties(ComputationalGraphAssembler assembler, TrainModel trainer) throws IOException {
        ModelPropertiesHelper mpHelper = new ModelPropertiesHelper();
        appendProperties(assembler, mpHelper);
        mpHelper.addProperties(getReaderProperties(trainer.args().trainingSets.get(0)));
        mpHelper.writeProperties(directory);
    }


    protected static int numLabels(INDArray labels) {
        FloatSet set = new FloatArraySet();
        for (int i = 0; i < labels.size(0); i++) {
            set.add(labels.getFloat(i));
        }
        return set.size();
    }


    public void appendProperties(ComputationalGraphAssembler assembler, ModelPropertiesHelper helper) {
        // give a chance to the assembler to save information to the model properties to describe the architecture
        // used for training:
        assembler.saveProperties(helper);

        //save the rest of the arguments:
        helper.setFeatureCalculator(featureMapper);
        helper.setLearningRate(args().learningRate);
        helper.setDropoutRate(args().dropoutRate);
        helper.setMiniBatchSize(args().miniBatchSize);
        // mpHelper.setBestScore(bestScore);
        helper.setNumEpochs(args().maxEpochs);
        helper.setNumTrainingSets(args().trainingSets.size());
        helper.setTime(time);
        helper.setSeed(args().seed);

        helper.setEarlyStopCriterion(args().stopWhenEpochsWithoutImprovement);
        helper.setRegularization(args().regularizationRate);
        helper.setPrecision(precision);
    }

    ParameterPrecision precision=ParameterPrecision.FP32;

    private static Properties getReaderProperties(String trainingSet) throws IOException {
        SequenceBaseInformationReader reader = new SequenceBaseInformationReader(trainingSet);
        final Properties properties = reader.getProperties();
        reader.close();
        return properties;
    }
}
