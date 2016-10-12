package org.campagnelab.dl.varanalysis.learning;

import com.beust.jcommander.JCommander;
import com.beust.jcommander.ParameterException;
import it.unimi.dsi.fastutil.floats.FloatArraySet;
import it.unimi.dsi.fastutil.floats.FloatSet;
import it.unimi.dsi.fastutil.objects.ObjectArrayList;
import org.apache.commons.io.FileUtils;
import org.apache.commons.io.FilenameUtils;
import org.campagnelab.dl.model.utils.mappers.FeatureMapper;
import org.campagnelab.dl.model.utils.mappers.LabelMapper;
import org.campagnelab.dl.model.utils.mappers.SimpleFeatureCalculator;
import org.campagnelab.dl.model.utils.models.ModelLoader;
import org.campagnelab.dl.varanalysis.learning.architecture.NeuralNetAssembler;
import org.campagnelab.dl.varanalysis.learning.architecture.SixDenseLayers;
import org.campagnelab.dl.varanalysis.learning.architecture.SixDenseLayersNarrower2;
import org.campagnelab.dl.varanalysis.learning.iterators.BaseInformationConcatIterator;
import org.campagnelab.dl.varanalysis.learning.iterators.BaseInformationIterator;
import org.campagnelab.dl.varanalysis.learning.iterators.FirstNIterator;
import org.campagnelab.dl.varanalysis.learning.models.ModelPropertiesHelper;
import org.campagnelab.dl.varanalysis.learning.models.PerformanceLogger;
import org.deeplearning4j.earlystopping.EarlyStoppingResult;
import org.deeplearning4j.earlystopping.saver.LocalFileModelSaver;
import org.deeplearning4j.nn.api.Layer;
import org.deeplearning4j.nn.api.Updater;
import org.deeplearning4j.nn.conf.LearningRatePolicy;
import org.deeplearning4j.nn.conf.MultiLayerConfiguration;
import org.deeplearning4j.nn.multilayer.MultiLayerNetwork;
import org.deeplearning4j.nn.weights.WeightInit;
import org.nd4j.jita.conf.CudaEnvironment;
import org.nd4j.linalg.api.ndarray.INDArray;
import org.nd4j.linalg.dataset.api.iterator.DataSetIterator;
import org.nd4j.linalg.factory.Nd4j;
import org.nd4j.linalg.lossfunctions.LossFunctions;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.io.*;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.Date;
import java.util.List;

/**
 * Abatract class to facilitate variations of training protocols.
 * Created by fac2003 on 7/12/16.
 */
public abstract class SomaticTrainer {
    static private Logger LOG = LoggerFactory.getLogger(TrainSomaticModel.class);
    protected TrainingArguments arguments;

    protected static TrainingArguments parseArguments(String[] args, String commandName) {

        TrainingArguments arguments = new TrainingArguments();
        JCommander commander = new JCommander(arguments);
        commander.setProgramName(commandName);
        try {
            commander.parse(args);
        } catch (ParameterException e) {

            commander.usage();
            throw e;

        }
        if (arguments.previousModelPath != null) {
            System.out.println(String.format("Resuming training with %s model parameters from %s %n", arguments.previousModelName, arguments.previousModelPath));
        }
        return arguments;
    }


    protected double dropoutRate = 0.5;
    protected LabelMapper labelMapper = new SimpleFeatureCalculator();
    protected FeatureMapper featureCalculator;
    protected String directory;
    protected long time;
    protected int numHiddenNodes;
    protected String attempt;
    protected double bestScore;
    protected int numTrainingFiles;
    protected MultiLayerNetwork net;
    protected LossFunctions.LossFunction lossFunction;
    protected PerformanceLogger performanceLogger;

    public SomaticTrainer(TrainingArguments arguments) {
        this.arguments = arguments;
    }

    public void execute(FeatureMapper featureCalculator, String trainingDataset[], int miniBatchSize) throws IOException {

        this.featureCalculator = featureCalculator;
        this.numTrainingFiles = trainingDataset.length;

        String path = "";

        time = new Date().getTime();
        System.out.println("time: " + time);
        System.out.println("epochs: " + arguments.maxEpochs);
        System.out.println(featureCalculator.getClass().getTypeName());
        directory = "models/" + Long.toString(time);
        attempt = "batch=" + miniBatchSize + "-learningRate=" + arguments.learningRate + "-time=" + time;
        int generateSamplesEveryNMinibatches = 10;
        FileUtils.forceMkdir(new File(directory));

        // Assemble the training iterator:
        labelMapper = new SimpleFeatureCalculator();
        List<BaseInformationIterator> trainIterList = new ObjectArrayList<>(trainingDataset.length);
        for (int i = 0; i < trainingDataset.length; i++) {
            trainIterList.add(new BaseInformationIterator(trainingDataset[i], miniBatchSize,
                    featureCalculator, labelMapper));
        }
        DataSetIterator async = new BaseInformationConcatIterator(trainIterList, miniBatchSize, featureCalculator, labelMapper);
        if (arguments.numTraining != Integer.MAX_VALUE) {
            async = new FirstNIterator(async, arguments.numTraining);
        }
        async = decorateIterator(async);
        System.out.println("Estimating scaling parameters:");
        //Load the training data:
        int numInputs = async.inputColumns();
        int numOutputs = async.totalOutcomes();
        numHiddenNodes = numInputs * 5;
        NeuralNetAssembler assembler = new SixDenseLayers();
        assembler.setSeed(arguments.seed);
        assembler.setLearningRate(arguments.learningRate);
        assembler.setNumHiddenNodes(numHiddenNodes);
        assembler.setNumInputs(numInputs);
        assembler.setNumOutputs(numOutputs);

        lossFunction = LossFunctions.LossFunction.MCXENT;

        assembler.setLossFunction(lossFunction);
        if (arguments.regularizationRate != Double.NaN) {
            assembler.setRegularization(true);
            assembler.setRegularizationRate(arguments.regularizationRate);
        }
        //   assembler.setDropoutRate(dropoutRate);


        //changed from XAVIER in iteration 14
        assembler.setWeightInitialization(WeightInit.RELU);
        assembler.setLearningRatePolicy(LearningRatePolicy.Score);
        MultiLayerConfiguration conf = assembler.createNetwork();
        net = new MultiLayerNetwork(conf);
        net.init();
        if (arguments.previousModelPath != null) {
            // Load the parameters of a previously trained model and set them on the new model to continue
            // training where we left it off. Note that models must have the same architecture or setting
            // parameters will fail.
            ModelLoader loader = new ModelLoader(arguments.previousModelPath);
            MultiLayerNetwork savedNetwork = loader.loadModel(arguments.previousModelName);
            if (savedNetwork == null || savedNetwork.getUpdater()==null || savedNetwork.params()==null) {
                System.err.println("Unable to load model or updater from " + arguments.previousModelPath);
            } else {
                net.setUpdater(savedNetwork.getUpdater());
                net.setParams(savedNetwork.params());

            }
        }

        //Print the  number of parameters in the network (and for each layer)
        Layer[] layers = net.getLayers();
        int totalNumParams = 0;
        for (int i = 0; i < layers.length; i++) {
            int nParams = layers[i].numParams();
            System.out.println("Number of parameters in layer " + i + ": " + nParams);
            totalNumParams += nParams;
        }
        System.out.println("Total number of network parameters: " + totalNumParams);

        writeProperties(this);
        performanceLogger = new PerformanceLogger(directory);
        EarlyStoppingResult<MultiLayerNetwork> result = train(conf, async);


        //Print out the results:
        System.out.println("Termination reason: " + result.getTerminationReason());
        System.out.println("Termination details: " + result.getTerminationDetails());
        System.out.println("Total epochs: " + result.getTotalEpochs());
        System.out.println("Best epoch number: " + result.getBestModelEpoch());
        System.out.println("Score at best epoch: " + result.getBestModelScore());

        writeProperties(this);
        writeBestScoreFile();
        System.out.println("Model completed, saved at time: " + attempt);
        performanceLogger.write();
    }

    protected DataSetIterator decorateIterator(DataSetIterator iterator) {
        return iterator;
    }

    protected void writeBestScoreFile() throws IOException {

        FileWriter scoreWriter = new FileWriter(directory + "/bestScore");
        scoreWriter.append(Double.toString(bestScore));
        scoreWriter.close();
    }

    protected void writeProperties(SomaticTrainer trainer) throws IOException {
        ModelPropertiesHelper mpHelper = new ModelPropertiesHelper(this);
        mpHelper.writeProperties(directory);
    }

    protected abstract EarlyStoppingResult<MultiLayerNetwork> train(MultiLayerConfiguration conf, DataSetIterator async)
            throws IOException;

    protected static void saveModel(LocalFileModelSaver saver, String directory, String prefix, MultiLayerNetwork net) throws IOException {
        FilenameUtils.concat(directory, prefix + "ModelConf.json");
        String confOut = FilenameUtils.concat(directory, prefix + "ModelConf.json");
        String paramOut = FilenameUtils.concat(directory, prefix + "ModelParams.bin");
        String updaterOut = FilenameUtils.concat(directory, prefix + "ModelUpdater.bin");
        save(net, confOut, paramOut, updaterOut);
    }

    protected static void save(MultiLayerNetwork net, String confOut, String paramOut, String updaterOut) throws IOException {
        String confJSON = net.getLayerWiseConfigurations().toJson();
        INDArray params = net.params();
        Updater updater = net.getUpdater();

        FileUtils.writeStringToFile(new File(confOut), confJSON, "UTF-8");
        try (DataOutputStream dos = new DataOutputStream(new BufferedOutputStream(Files.newOutputStream(Paths.get(paramOut))))) {
            Nd4j.write(params, dos);
        }

        try (ObjectOutputStream oos = new ObjectOutputStream(new BufferedOutputStream(new FileOutputStream(new File(updaterOut))))) {
            oos.writeObject(updater);
        }
    }

    protected static int numLabels(INDArray labels) {
        FloatSet set = new FloatArraySet();
        for (int i = 0; i < labels.size(0); i++) {
            set.add(labels.getFloat(i));
        }
        return set.size();
    }

    public void appendProperties(ModelPropertiesHelper helper) {
        helper.setFeatureCalculator(featureCalculator);
        helper.setLearningRate(arguments.learningRate);
        helper.setNumHiddenNodes(numHiddenNodes);
        helper.setMiniBatchSize(arguments.miniBatchSize);
        // mpHelper.setBestScore(bestScore);
        helper.setNumEpochs(arguments.maxEpochs);
        helper.setNumTrainingSets(numTrainingFiles);
        helper.setTime(time);
        helper.setSeed(arguments.seed);
        helper.setLossFunction(lossFunction.name());
        helper.setEarlyStopCriterion(arguments.stopWhenEpochsWithoutImprovement);
        helper.setRegularization(arguments.regularizationRate);
    }
}
