package org.campagnelab.dl.varanalysis.learning;

import com.beust.jcommander.JCommander;
import com.beust.jcommander.ParameterException;
import it.unimi.dsi.fastutil.floats.FloatArraySet;
import it.unimi.dsi.fastutil.floats.FloatSet;
import it.unimi.dsi.fastutil.objects.ObjectArrayList;
import org.apache.commons.io.FileUtils;
import org.apache.commons.io.FilenameUtils;
import org.campagnelab.dl.model.utils.ConfigurableFeatureMapper;
import org.campagnelab.dl.model.utils.mappers.FeatureMapper;
import org.campagnelab.dl.model.utils.mappers.LabelMapper;
import org.campagnelab.dl.model.utils.mappers.SimpleFeatureCalculator;
import org.campagnelab.dl.model.utils.models.ModelLoader;
import org.campagnelab.dl.varanalysis.learning.architecture.NeuralNetAssembler;
import org.campagnelab.dl.varanalysis.learning.iterators.BaseInformationConcatIterator;
import org.campagnelab.dl.varanalysis.learning.iterators.BaseInformationIterator;
import org.campagnelab.dl.varanalysis.learning.iterators.FirstNIterator;
import org.campagnelab.dl.varanalysis.learning.iterators.NamedDataSetIterator;
import org.campagnelab.dl.varanalysis.learning.models.ModelPropertiesHelper;
import org.campagnelab.dl.varanalysis.learning.models.PerformanceLogger;
import org.campagnelab.dl.varanalysis.tools.AbstractTool;
import org.campagnelab.goby.baseinfo.SequenceBaseInformationReader;
import org.deeplearning4j.earlystopping.EarlyStoppingResult;
import org.deeplearning4j.earlystopping.saver.LocalFileModelSaver;
import org.deeplearning4j.nn.api.Layer;
import org.deeplearning4j.nn.api.Updater;
import org.deeplearning4j.nn.conf.LearningRatePolicy;
import org.deeplearning4j.nn.conf.MultiLayerConfiguration;
import org.deeplearning4j.nn.multilayer.MultiLayerNetwork;
import org.deeplearning4j.nn.weights.WeightInit;
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
import java.util.Properties;

/**
 * Abstract class to facilitate variations of training protocols.
 * Created by fac2003 on 7/12/16.
 */
public abstract class SomaticTrainer extends AbstractTool<TrainingArguments> {
    static private Logger LOG = LoggerFactory.getLogger(TrainSomaticModel.class);

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
        return arguments;
    }

    @Override
    public void execute() {
        if (args().getTrainingSets().length == 0) {
            System.err.println("You must provide training datasets.");
        }
        FeatureMapper featureMapper = null;
        try {
            featureMapper = configureFeatureMapper(args().featureMapperClassname, args().isTrio, args().getTrainingSets());
            execute(featureMapper, args().getTrainingSets(), args().miniBatchSize);
        } catch (IOException e) {
            System.err.println("An exception occured. Details may be provided below");
            e.printStackTrace();
        }

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


    public void execute(FeatureMapper featureCalculator, String trainingDataset[], int miniBatchSize) throws IOException {
        if (args().previousModelPath != null) {
            System.out.println(String.format("Resuming training with %s model parameters from %s %n", args().previousModelName, args().previousModelPath));
        }
        this.featureCalculator = featureCalculator;
        this.numTrainingFiles = trainingDataset.length;

        String path = "";

        time = new Date().getTime();
        System.out.println("time: " + time);
        System.out.println("epochs: " + args().maxEpochs);
        System.out.println(featureCalculator.getClass().getTypeName());
        directory = "models/" + Long.toString(time);
        attempt = "batch=" + miniBatchSize + "-learningRate=" + args().learningRate + "-time=" + time;
        int generateSamplesEveryNMinibatches = 10;
        FileUtils.forceMkdir(new File(directory));

        // Assemble the training iterator:
        labelMapper = new SimpleFeatureCalculator();
        List<BaseInformationIterator> trainIterList = new ObjectArrayList<>(trainingDataset.length);
        for (int i = 0; i < trainingDataset.length; i++) {
            trainIterList.add(new BaseInformationIterator(trainingDataset[i], miniBatchSize,
                    featureCalculator, labelMapper));
        }
        NamedDataSetIterator async = new BaseInformationConcatIterator(trainIterList, miniBatchSize, featureCalculator, labelMapper);
        if (args().numTraining != Integer.MAX_VALUE) {
            async = new FirstNIterator(async, args().numTraining);
        }
        async = decorateIterator(async);
        System.out.println("Estimating scaling parameters:");
        //Load the training data:
        int numInputs = async.inputColumns();
        int numOutputs = async.totalOutcomes();
        numHiddenNodes = numInputs * 5;
        NeuralNetAssembler assembler = getNeuralNetAssembler();
        assembler.setSeed(args().seed);
        assembler.setLearningRate(args().learningRate);
        assembler.setDropoutRate(args().dropoutRate);
        assembler.setNumHiddenNodes(numHiddenNodes);
        assembler.setNumInputs(numInputs);
        assembler.setNumOutputs(numOutputs);

        lossFunction = LossFunctions.LossFunction.MCXENT;

        assembler.setLossFunction(lossFunction);
        if (args().regularizationRate != Double.NaN) {
            assembler.setRegularization(true);
            assembler.setRegularizationRate(args().regularizationRate);
        }
        //   assembler.setDropoutRate(dropoutRate);


        //changed from XAVIER in iteration 14
        assembler.setWeightInitialization(WeightInit.RELU);
        assembler.setLearningRatePolicy(LearningRatePolicy.Score);
        MultiLayerConfiguration conf = assembler.createNetwork();
        net = new MultiLayerNetwork(conf);
        net.init();
        if (args().previousModelPath != null) {
            // Load the parameters of a previously trained model and set them on the new model to continue
            // training where we left it off. Note that models must have the same architecture or setting
            // parameters will fail.
            ModelLoader loader = new ModelLoader(args().previousModelPath);
            MultiLayerNetwork savedNetwork = loader.loadModel(args().previousModelName);
            if (savedNetwork == null || savedNetwork.getUpdater() == null || savedNetwork.params() == null) {
                System.err.println("Unable to load model or updater from " + args().previousModelPath);
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
        System.out.println("Score at best epoch: " + performanceLogger.getBestScore());
        System.out.println("AUC at best epoch: " + performanceLogger.getBestAUC());

        writeProperties(this);
        writeBestScoreFile();
        System.out.println("Model completed, saved at time: " + attempt);
        performanceLogger.write();
    }

    private NeuralNetAssembler getNeuralNetAssembler() {
        try {
            return (NeuralNetAssembler) Class.forName(args().architectureClassname).newInstance();
        } catch (Exception e) {
            System.err.println("Unable to instantiate net architecture " + args().architectureClassname);
            System.exit(1);
        }
        return null;
    }


    protected NamedDataSetIterator decorateIterator(NamedDataSetIterator iterator) {
        return iterator;
    }

    protected void writeBestScoreFile() throws IOException {

        FileWriter scoreWriter = new FileWriter(directory + "/bestScore");
        scoreWriter.append(Double.toString(bestScore));
        scoreWriter.close();
    }

    protected void writeProperties(SomaticTrainer trainer) throws IOException {
        ModelPropertiesHelper mpHelper = new ModelPropertiesHelper(this);
        mpHelper.addProperties(getReaderProperties(trainer.args().trainingSets.get(0)));
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

    protected Precision precision = Precision.FP32;

    public enum Precision {
        FP16,
        FP32
    }

    public void appendProperties(ModelPropertiesHelper helper) {
        helper.setFeatureCalculator(featureCalculator);
        helper.setLearningRate(args().learningRate);
        helper.setDropoutRate(args().dropoutRate);
        helper.setNumHiddenNodes(numHiddenNodes);
        helper.setMiniBatchSize(args().miniBatchSize);
        // mpHelper.setBestScore(bestScore);
        helper.setNumEpochs(args().maxEpochs);
        helper.setNumTrainingSets(numTrainingFiles);
        helper.setTime(time);
        helper.setSeed(args().seed);
        helper.setLossFunction(lossFunction.name());
        helper.setEarlyStopCriterion(args().stopWhenEpochsWithoutImprovement);
        helper.setRegularization(args().regularizationRate);
        helper.setPrecision(precision);
    }


    public static FeatureMapper configureFeatureMapper(String featureMapperClassname, boolean isTrio, String[] trainingSets) throws IOException {


        try {
            Class clazz = Class.forName(featureMapperClassname + (isTrio ? "Trio" : ""));
            final FeatureMapper featureMapper = (FeatureMapper) clazz.newInstance();
            if (featureMapper instanceof ConfigurableFeatureMapper) {
                ConfigurableFeatureMapper cmapper = (ConfigurableFeatureMapper) featureMapper;
                if (trainingSets.length > 1) {
                    LOG.warn("sbip properties are only read from the first training set. Concat the files before training if you need to use properties across all inputs.");
                }
                final Properties properties = getReaderProperties(trainingSets[0]);
                cmapper.configure(properties);
            }
            return featureMapper;
        } catch (IllegalAccessException e) {
            e.printStackTrace();
        } catch (InstantiationException e) {
            e.printStackTrace();
        } catch (ClassNotFoundException e) {
            e.printStackTrace();
        }
        return null;
    }

    private static Properties getReaderProperties(String trainingSet) throws IOException {
        SequenceBaseInformationReader reader = new SequenceBaseInformationReader(trainingSet);
        final Properties properties = reader.getProperties();
        reader.close();
        return properties;
    }
}
