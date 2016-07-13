package org.campagnelab.dl.varanalysis.learning;

import it.unimi.dsi.fastutil.objects.ObjectArrayList;
import org.apache.commons.io.FileUtils;
import org.apache.commons.io.FilenameUtils;
import org.campagnelab.dl.model.utils.mappers.FeatureMapper;
import org.campagnelab.dl.model.utils.mappers.LabelMapper;
import org.campagnelab.dl.model.utils.mappers.SimpleFeatureCalculator;
import org.campagnelab.dl.varanalysis.learning.architecture.NeuralNetAssembler;
import org.campagnelab.dl.varanalysis.learning.architecture.SixDenseLayersNarrower2;
import org.campagnelab.dl.varanalysis.learning.iterators.BaseInformationConcatIterator;
import org.campagnelab.dl.varanalysis.learning.iterators.BaseInformationIterator;
import org.deeplearning4j.datasets.iterator.AsyncDataSetIterator;
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
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.io.*;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.Date;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

/**
 * Abatract class to facilitate variations of training protocols.
 * Created by fac2003 on 7/12/16.
 */
public abstract class SomaticTrainer {
    static private Logger LOG = LoggerFactory.getLogger(TrainSomaticModel.class);
    protected int seed = 123;
    protected double learningRate = 0.1;
    protected int miniBatchSize = 32;
    protected int numEpochs = 15;
    protected int earlyStopCondition = 3;
    protected double dropoutRate = 0.5;
    protected LabelMapper labelMapper;
    protected FeatureMapper featureCalculator;
    protected String directory;
    protected long time;
    protected int numHiddenNodes;
    protected String attempt;
    protected double bestScore;
    protected int numTrainingFiles;
    protected MultiLayerNetwork net;
    public void execute(FeatureMapper featureCalculator, String trainingDataset[], int miniBatchSize) throws IOException {
        this.featureCalculator = featureCalculator;
this.miniBatchSize=miniBatchSize;

        String path = "";

        time = new Date().getTime();
        System.out.println("time: " + time);
        System.out.println("epochs: " + numEpochs);
        System.out.println(featureCalculator.getClass().getTypeName());
        directory = "models/" + Long.toString(time);
        attempt = "batch=" + miniBatchSize + "-learningRate=" + learningRate + "-time=" + time;
        int generateSamplesEveryNMinibatches = 10;
        FileUtils.forceMkdir(new File(directory));

        // Assemble the training iterator:
        labelMapper = new SimpleFeatureCalculator();
        List<BaseInformationIterator> trainIterList = new ObjectArrayList<>(trainingDataset.length);
        for (int i = 0; i < trainingDataset.length; i++) {
            trainIterList.add(new BaseInformationIterator(trainingDataset[i], miniBatchSize,
                    featureCalculator, labelMapper));
        }
        final AsyncDataSetIterator async = new AsyncDataSetIterator(new BaseInformationConcatIterator(trainIterList, miniBatchSize*2, featureCalculator, labelMapper));

        System.out.println("Estimating scaling parameters:");
        //Load the training data:
        int numInputs = async.inputColumns();
        int numOutputs = async.totalOutcomes();
        numHiddenNodes = numInputs * 5;
        NeuralNetAssembler assembler = new SixDenseLayersNarrower2();
        assembler.setSeed(seed);
        assembler.setLearningRate(learningRate);
        assembler.setNumHiddenNodes(numHiddenNodes);
        assembler.setNumInputs(numInputs);
        assembler.setNumOutputs(numOutputs);
        assembler.setRegularization(false);
        // assembler.setRegularizationRate(1e-6);
        // assembler.setDropoutRate(dropoutRate);


        //changed from XAVIER in iteration 14
        assembler.setWeightInitialization(WeightInit.RELU);
        assembler.setLearningRatePolicy(LearningRatePolicy.Score);
        MultiLayerConfiguration conf = assembler.createNetwork();
        net= new MultiLayerNetwork(conf);
        net.init();
        //Print the  number of parameters in the network (and for each layer)
        Layer[] layers = net.getLayers();
        int totalNumParams = 0;
        for (int i = 0; i < layers.length; i++) {
            int nParams = layers[i].numParams();
            System.out.println("Number of parameters in layer " + i + ": " + nParams);
            totalNumParams += nParams;
        }
        System.out.println("Total number of network parameters: " + totalNumParams);

        writeProperties(featureCalculator);

        EarlyStoppingResult<MultiLayerNetwork> result = train(conf, async);


        //Print out the results:
        System.out.println("Termination reason: " + result.getTerminationReason());
        System.out.println("Termination details: " + result.getTerminationDetails());
        System.out.println("Total epochs: " + result.getTotalEpochs());
        System.out.println("Best epoch number: " + result.getBestModelEpoch());
        System.out.println("Score at best epoch: " + result.getBestModelScore());

        writeProperties(featureCalculator);
        writeBestScoreFile();
        System.out.println("Model completed, saved at time: " + attempt);

    }

    protected void writeBestScoreFile() throws IOException {

        FileWriter scoreWriter = new FileWriter(directory + "/bestScore");
        scoreWriter.append(Double.toString(bestScore));
        scoreWriter.close();
    }

    protected void writeProperties(FeatureMapper featureCalculator) throws IOException {
        ModelPropertiesHelper mpHelper = new ModelPropertiesHelper();
        mpHelper.setFeatureCalculator(featureCalculator);
        mpHelper.setLearningRate(learningRate);
        mpHelper.setNumHiddenNodes(numHiddenNodes);
        mpHelper.setMiniBatchSize(miniBatchSize);
        // mpHelper.setBestScore(bestScore);
        mpHelper.setNumEpochs(numEpochs);
        mpHelper.setNumTrainingSets(numTrainingFiles);
        mpHelper.setTime(time);
        mpHelper.setSeed(seed);
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
        Set<Float> set = new HashSet<>();
        //for (labels.size(1));
        return 2;
    }

}
