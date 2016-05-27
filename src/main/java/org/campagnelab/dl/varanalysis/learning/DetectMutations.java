package org.campagnelab.dl.varanalysis.learning;

import it.unimi.dsi.logging.ProgressLogger;
import org.campagnelab.dl.varanalysis.learning.iterators.*;
import org.campagnelab.dl.varanalysis.learning.mappers.FeatureMapper;
import org.campagnelab.dl.varanalysis.learning.mappers.FeatureMapperV2;
import org.campagnelab.dl.varanalysis.learning.mappers.LabelMapper;
import org.campagnelab.dl.varanalysis.learning.mappers.PositiveControlFeatureMapper;
import org.deeplearning4j.earlystopping.saver.LocalFileModelSaver;
import org.deeplearning4j.nn.api.Layer;
import org.deeplearning4j.nn.api.OptimizationAlgorithm;
import org.deeplearning4j.nn.conf.MultiLayerConfiguration;
import org.deeplearning4j.nn.conf.NeuralNetConfiguration;
import org.deeplearning4j.nn.conf.Updater;
import org.deeplearning4j.nn.conf.layers.DenseLayer;
import org.deeplearning4j.nn.conf.layers.OutputLayer;
import org.deeplearning4j.nn.multilayer.MultiLayerNetwork;
import org.deeplearning4j.nn.weights.WeightInit;
import org.deeplearning4j.optimize.api.IterationListener;
import org.deeplearning4j.optimize.listeners.ScoreIterationListener;
import org.nd4j.linalg.api.ndarray.INDArray;
import org.nd4j.linalg.dataset.DataSet;
import org.nd4j.linalg.lossfunctions.LossFunctions;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.io.*;
import java.util.Collections;
import java.util.Date;
import java.util.HashSet;
import java.util.Set;

/**
 * Train a neural network to predict mutations.
 * <p>
 * Created by fac2003 on 5/21/16.
 *
 * @author Fabien Campagne
 */
public class DetectMutations {
    static private Logger LOG = LoggerFactory.getLogger(DetectMutations.class);


    public static void main(String[] args) throws IOException {
        int seed = 123;
        double learningRate = 0.05;
        int miniBatchSize = 100;
        int numEpochs = 3;
        String attempt = "batch=" + miniBatchSize + "learningRate=" + learningRate + "-time:" + new Date().getTime();
        int generateSamplesEveryNMinibatches = 10;

        //Load the training data:
        final FeatureMapper featureCalculator = new FeatureMapperV2();//new PositiveControlFeatureMapper();//
        final LabelMapper labelMapper = new SimpleFeatureCalculator();
        BaseInformationIterator trainIter = new BaseInformationIterator("sample_data/protobuf/genotypes_proto_mutated_randomized.parquet",
                miniBatchSize, featureCalculator, labelMapper);

        int numInputs = trainIter.inputColumns();
        int numOutputs = trainIter.totalOutcomes();
        int numHiddenNodes = 500;

        MultiLayerConfiguration conf = new NeuralNetConfiguration.Builder()
                .seed(seed)
                .iterations(10)
                .optimizationAlgo(OptimizationAlgorithm.STOCHASTIC_GRADIENT_DESCENT)
                .learningRate(learningRate).regularization(true).l2(0.00002)
                .updater(Updater.ADAGRAD)
                .list()
                .layer(0, new DenseLayer.Builder().nIn(numInputs).nOut(numHiddenNodes)
                        .weightInit(WeightInit.XAVIER)
                        .activation("relu")
                        .build())
                .layer(1, new DenseLayer.Builder().nIn(numHiddenNodes).nOut(numHiddenNodes)
                        .weightInit(WeightInit.XAVIER)
                        .activation("relu")
                        .build())
                .layer(2, new OutputLayer.Builder(LossFunctions.LossFunction.NEGATIVELOGLIKELIHOOD)
                        .weightInit(WeightInit.XAVIER)
                        .activation("softmax").weightInit(WeightInit.XAVIER)
                        .nIn(numHiddenNodes).nOut(numOutputs).build())
                .pretrain(false).backprop(true).build();


        MultiLayerNetwork net = new MultiLayerNetwork(conf);
        net.init();
    //    net.setListeners(/*new HistogramIterationListener(1), */new ScoreIterationListener(1));
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

        for (int i = 0; i < numEpochs; i++) {
            ProgressLogger pg = new ProgressLogger(LOG);
            pg.itemsName = "mini-batch";

            pg.expectedUpdates = trainIter.totalExamples() / miniBatchSize; // one iteration processes miniBatchIterator elements.
            pg.start();
            while (trainIter.hasNext()) {
                // trigger bugs early:
//                printSample(miniBatchSize, exampleLength, nSamplesToGenerate, nCharactersToSample, generationInitialization, rng, attempt, iter, net, miniBatchNumber);

                DataSet ds = trainIter.next(miniBatchSize);
                if (numLabels(ds.getLabels()) != 2) {
                    System.out.println("There should be two labels in the miniBatch");
                }

                net.fit(ds);

                INDArray predictedLabels = net.output(ds.getFeatures(), false);


                pg.update();

                if (net.score() < bestScore) {
                    bestScore = net.score();
                    saver.saveBestModel(net, net.score());
                    System.out.println("Saving best score model.. score=" + bestScore);
                }

            }
            pg.stop();
            pgEpoch.update();
            trainIter.reset();    //Reset iterator for another epoch
        }
        pgEpoch.stop();
        System.out.println("Saving last model with score=" + net.score());
        saver.saveLatestModel(net, net.score());

    }

    private static int numLabels(INDArray labels) {
        Set<Float> set = new HashSet<>();
        //for (labels.size(1));
        return 2;
    }

}
