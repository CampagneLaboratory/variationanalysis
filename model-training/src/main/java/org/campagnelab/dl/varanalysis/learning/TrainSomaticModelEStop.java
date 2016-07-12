package org.campagnelab.dl.varanalysis.learning;

import it.unimi.dsi.fastutil.objects.ObjectArrayList;
import org.apache.commons.io.FileUtils;
import org.apache.commons.io.FilenameUtils;
import org.campagnelab.dl.model.utils.mappers.*;
import org.campagnelab.dl.varanalysis.learning.architecture.*;
import org.campagnelab.dl.varanalysis.learning.iterators.*;
import org.campagnelab.dl.varanalysis.storage.RecordWriter;
import org.deeplearning4j.datasets.iterator.AsyncDataSetIterator;
import org.deeplearning4j.earlystopping.EarlyStoppingConfiguration;
import org.deeplearning4j.earlystopping.EarlyStoppingResult;
import org.deeplearning4j.earlystopping.saver.LocalFileModelSaver;
import org.deeplearning4j.earlystopping.scorecalc.DataSetLossCalculator;
import org.deeplearning4j.earlystopping.termination.MaxEpochsTerminationCondition;
import org.deeplearning4j.earlystopping.termination.ScoreImprovementEpochTerminationCondition;
import org.deeplearning4j.earlystopping.trainer.EarlyStoppingTrainer;
import org.deeplearning4j.nn.api.Updater;
import org.deeplearning4j.nn.conf.LearningRatePolicy;
import org.deeplearning4j.nn.conf.MultiLayerConfiguration;
import org.deeplearning4j.nn.multilayer.MultiLayerNetwork;
import org.deeplearning4j.nn.weights.WeightInit;
import org.nd4j.linalg.api.ndarray.INDArray;
import org.nd4j.linalg.factory.Nd4j;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.io.*;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.*;

/**
 * Train a neural network to predict mutations.
 * <p>
 * Created by fac2003 on 5/21/16.
 *
 * @author Fabien Campagne
 */
public class TrainSomaticModelEStop {
    public static final int MIN_ITERATION_BETWEEN_BEST_MODEL = 1000;
    static private Logger LOG = LoggerFactory.getLogger(TrainSomaticModel.class);


    public static void main(String[] args) throws IOException {
        final FeatureMapper featureCalculator = new FeatureMapperV13();
        if (args.length < 2) {
            System.err.println("usage: DetectMutations <input-validation-file> <input-training-directory>");
        }
        //VALIDATION FILE IS FIRST ARGUMENT
        String valFile = RecordWriter.addParqExtension(args[0]);
        int i = 0;
        String path = "";
        int seed = 123;
        double learningRate = 0.1;
        int miniBatchSize = 100;
        int numEpochs = 1;
        int earlyStopCondition = 3;
        double dropoutRate = 0.5;
        long time = new Date().getTime();
        System.out.println("time: " + time);
        System.out.println("epochs: " + numEpochs);
        System.out.println(featureCalculator.getClass().getTypeName());
        String directory = "models/" + Long.toString(time);
        String attempt = "batch=" + miniBatchSize + "-learningRate=" + learningRate + "-time=" + time;
        int generateSamplesEveryNMinibatches = 10;
        FileUtils.forceMkdir(new File(directory));



        System.out.println("Estimating scaling parameters:");
        final LabelMapper labelMapper = new SimpleFeatureCalculator();
        List<BaseInformationIterator> trainIterList = new ObjectArrayList<>(args.length-1);
        for (i = 1; i < args.length; i++){
            trainIterList.add(new BaseInformationIterator(RecordWriter.addParqExtension(args[i]), miniBatchSize,
                    featureCalculator, labelMapper));
        }
        final BaseInformationConcatIterator trainIter = new BaseInformationConcatIterator(trainIterList, miniBatchSize,
                featureCalculator, labelMapper);
        final AsyncDataSetIterator async = new AsyncDataSetIterator(trainIter);
        //Load the training data:
        int numInputs = async.inputColumns();
        int numOutputs = async.totalOutcomes();
        int numHiddenNodes = numInputs * 5;
        NeuralNetAssembler assembler = new SixDenseLayersNarrower2();
        assembler.setSeed(seed);
        assembler.setLearningRate(learningRate);
        assembler.setNumHiddenNodes(numHiddenNodes);
        assembler.setNumInputs(numInputs);
        assembler.setNumOutputs(numOutputs);
        assembler.setRegularization(false);
        // assembler.setRegularizationRate(1e-6);
        //   assembler.setDropoutRate(dropoutRate);


        //changed from XAVIER in iteration 14
        assembler.setWeightInitialization(WeightInit.RELU);
        assembler.setLearningRatePolicy(LearningRatePolicy.Score);
        MultiLayerConfiguration conf = assembler.createNetwork();

        //LocalFileModelSaver saver = new LocalFileModelSaver(attempt);

        EarlyStoppingConfiguration<MultiLayerNetwork> esConf = new EarlyStoppingConfiguration.Builder()
                .epochTerminationConditions(new MaxEpochsTerminationCondition(numEpochs),new ScoreImprovementEpochTerminationCondition(earlyStopCondition))
                .scoreCalculator(new DataSetLossCalculator(new BaseInformationIterator(valFile, miniBatchSize,
                featureCalculator, labelMapper), true))
                .evaluateEveryNEpochs(1)
                .modelSaver(new LocalFileModelSaver(directory))
                .build();

        EarlyStoppingTrainer trainer = new EarlyStoppingTrainer(esConf,conf,async);
        trainer.setListener(new EStatusListener());
        EarlyStoppingResult<MultiLayerNetwork> result = trainer.fit();


        //Print out the results:
        System.out.println("Termination reason: " + result.getTerminationReason());
        System.out.println("Termination details: " + result.getTerminationDetails());
        System.out.println("Total epochs: " + result.getTotalEpochs());
        System.out.println("Best epoch number: " + result.getBestModelEpoch());
        System.out.println("Score at best epoch: " + result.getBestModelScore());



        MultiLayerNetwork bestModel = result.getBestModel();
        double bestScore = bestModel.score();


        FileWriter scoreWriter = new FileWriter(directory + "/bestScore");
        scoreWriter.append(Double.toString(bestScore));
        scoreWriter.close();

        ModelPropertiesHelper mpHelper=new ModelPropertiesHelper();
        mpHelper.setFeatureCalculator(featureCalculator);
        mpHelper.setLearningRate(learningRate);
        mpHelper.setNumHiddenNodes(numHiddenNodes);
        mpHelper.setMiniBatchSize(miniBatchSize);
        mpHelper.setBestScore(bestScore);
        mpHelper.setNumEpochs(numEpochs);
        mpHelper.setNumTrainingSets(args.length);
        mpHelper.setTime(time);
        mpHelper.setSeed(seed);
        mpHelper.writeProperties(directory);


        System.out.println("Model completed, saved at time: " + attempt);

    }

    public static void saveModel(LocalFileModelSaver saver, String directory, String prefix, MultiLayerNetwork net) throws IOException {
        FilenameUtils.concat(directory, prefix + "ModelConf.json");
        String confOut = FilenameUtils.concat(directory, prefix + "ModelConf.json");
        String paramOut = FilenameUtils.concat(directory, prefix + "ModelParams.bin");
        String updaterOut = FilenameUtils.concat(directory, prefix + "ModelUpdater.bin");
        save(net, confOut, paramOut, updaterOut);
    }

    private static void save(MultiLayerNetwork net, String confOut, String paramOut, String updaterOut) throws IOException {
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

    private static int numLabels(INDArray labels) {
        Set<Float> set = new HashSet<>();
        //for (labels.size(1));
        return 2;
    }

}
