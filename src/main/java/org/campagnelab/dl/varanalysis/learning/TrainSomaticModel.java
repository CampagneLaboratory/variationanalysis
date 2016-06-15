package org.campagnelab.dl.varanalysis.learning;

import it.unimi.dsi.fastutil.objects.ObjectArrayList;
import it.unimi.dsi.logging.ProgressLogger;
import org.apache.commons.io.FileUtils;
import org.apache.commons.io.FilenameUtils;
import org.campagnelab.dl.varanalysis.learning.architecture.*;
import org.campagnelab.dl.varanalysis.learning.iterators.*;
import org.campagnelab.dl.varanalysis.learning.mappers.*;
import org.deeplearning4j.datasets.iterator.AsyncDataSetIterator;
import org.deeplearning4j.earlystopping.saver.LocalFileModelSaver;
import org.deeplearning4j.nn.api.Layer;
import org.deeplearning4j.nn.api.Updater;
import org.deeplearning4j.nn.conf.LearningRatePolicy;
import org.deeplearning4j.nn.conf.MultiLayerConfiguration;
import org.deeplearning4j.nn.multilayer.MultiLayerNetwork;
import org.deeplearning4j.nn.weights.WeightInit;
import org.nd4j.linalg.api.ndarray.INDArray;
import org.nd4j.linalg.dataset.DataSet;
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
 * Train a neural network to predict mutations.
 * <p>
 * Created by fac2003 on 5/21/16.
 *
 * @author Fabien Campagne
 */
public class TrainSomaticModel {
    public static final int MIN_ITERATION_BETWEEN_BEST_MODEL = 1000;
    static private Logger LOG = LoggerFactory.getLogger(TrainSomaticModel.class);


    public static void main(String[] args) throws IOException {
        final FeatureMapper featureCalculator = new FeatureMapperV9();
        if (args.length < 1) {
            System.err.println("usage: DetectMutations <input-training-file_1> <input-training-file_2> ... ");
        }


        int seed = 123;
        double learningRate = 0.1;
        int miniBatchSize = 100;
        int numEpochs = 50;
        double dropoutRate = 0.5;
        long time = new Date().getTime();
        System.out.println("time: " + time);
        System.out.println("epochs: " + numEpochs);
        System.out.println(featureCalculator.getClass().getTypeName());
        String attempt = "batch=" + miniBatchSize + "-learningRate=" + learningRate + "-time=" + time;
        int generateSamplesEveryNMinibatches = 10;
        FileUtils.forceMkdir(new File(attempt));

        System.out.println("Estimating scaling parameters:");
        final LabelMapper labelMapper = new SimpleFeatureCalculator();
        List<BaseInformationIterator> trainIterList = new ObjectArrayList<>(args.length);
        for (int i = 0; i < args.length; i++){
            trainIterList.add(new BaseInformationIterator(args[i], miniBatchSize,
                    featureCalculator, labelMapper));
        }
        final AsyncDataSetIterator async = new AsyncDataSetIterator(new BaseInformationConcatIterator(trainIterList,miniBatchSize,featureCalculator,labelMapper));
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


        MultiLayerNetwork net = new MultiLayerNetwork(conf);
        net.init();
        //      net.setListeners(new HistogramIterationListener(10) /*, new ScoreIterationListener(1) */);
        //Print the  number of parameters in the network (and for each layer)
        Layer[] layers = net.getLayers();
        int totalNumParams = 0;
        for (int i = 0; i < layers.length; i++) {
            int nParams = layers[i].numParams();
            System.out.println("Number of parameters in layer " + i + ": " + nParams);
            totalNumParams += nParams;
        }
        System.out.println("Total number of network parameters: " + totalNumParams);

        //Do training, and then generate and print samples from network
        int miniBatchNumber = 0;
        boolean init = true;
        ProgressLogger pgEpoch = new ProgressLogger(LOG);
        pgEpoch.itemsName = "epoch";
        pgEpoch.expectedUpdates = numEpochs;
        pgEpoch.start();
        double bestScore = Double.MAX_VALUE;
        LocalFileModelSaver saver = new LocalFileModelSaver(attempt);
        int iter = 0;
        for (int epoch = 0; epoch < numEpochs; epoch++) {
            ProgressLogger pg = new ProgressLogger(LOG);
            pg.itemsName = "mini-batch";

            pg.expectedUpdates = async.totalExamples() / miniBatchSize; // one iteration processes miniBatchIterator elements.
            pg.start();
            int lastIter=0;
            while (async.hasNext()) {
                // trigger bugs early:
//                printSample(miniBatchSize, exampleLength, nSamplesToGenerate, nCharactersToSample, generationInitialization, rng, attempt, iter, net, miniBatchNumber);

                DataSet ds = async.next();
                if (numLabels(ds.getLabels()) != 2) {
                    System.out.println("There should be two labels in the miniBatch");
                }
                // scale the features:
                //   scaler.transform(ds);

                // fit the net:
                net.fit(ds);
                INDArray predictedLabels = net.output(ds.getFeatures(), false);


                pg.update();

                double score = net.score();
                if (Double.isNaN(score)) {
                    //   System.out.println(net.params());
                    System.out.println("nan at " + iter);
                }
                iter++;
                if (score < bestScore * 0.95 && iter>(lastIter+ MIN_ITERATION_BETWEEN_BEST_MODEL)) {
                    bestScore = score;
                    saver.saveBestModel(net, score);
                    System.out.println("Saving best score model.. score=" + bestScore);
                    lastIter=iter;
                }

            }
            //save latest after the end of an epoch:
            saver.saveLatestModel(net, net.score());
            saveModel(saver, attempt, Integer.toString(epoch) + "-", net);

            pg.stop();
            pgEpoch.update();
            async.reset();    //Reset iterator for another epoch
        }
        pgEpoch.stop();
        System.out.println("Saving last model with score=" + net.score());
        saver.saveLatestModel(net, net.score());
        FileWriter scoreWriter = new FileWriter(attempt + "/bestScore");
        scoreWriter.append(Double.toString(bestScore));
        scoreWriter.close();

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
