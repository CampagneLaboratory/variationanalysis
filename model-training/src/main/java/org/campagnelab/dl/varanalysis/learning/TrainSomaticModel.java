package org.campagnelab.dl.varanalysis.learning;

import it.unimi.dsi.fastutil.objects.ObjectArrayList;
import it.unimi.dsi.logging.ProgressLogger;
import org.apache.commons.collections.map.HashedMap;
import org.apache.commons.io.FileUtils;
import org.apache.commons.io.FilenameUtils;
import org.campagnelab.dl.model.utils.mappers.*;
import org.campagnelab.dl.varanalysis.learning.architecture.*;
import org.campagnelab.dl.varanalysis.learning.iterators.*;
import org.campagnelab.dl.varanalysis.storage.RecordWriter;
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
import org.nd4j.linalg.dataset.DataSet;
import org.nd4j.linalg.dataset.api.iterator.DataSetIterator;
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
public class TrainSomaticModel extends SomaticTrainer {
    public static final int MIN_ITERATION_BETWEEN_BEST_MODEL = 1000;
    static private Logger LOG = LoggerFactory.getLogger(TrainSomaticModel.class);


    public static void main(String[] args) throws IOException {

        TrainSomaticModel trainer = new TrainSomaticModel();
        if (args.length < 1) {
            System.err.println("usage: DetectMutations <input-training-directory>");
        }
        trainer.execute(new FeatureMapperV15(), args);

    }


    @Override
    protected EarlyStoppingResult<MultiLayerNetwork> train(MultiLayerConfiguration conf, DataSetIterator async) throws IOException {
        //Do training, and then generate and print samples from network
        int miniBatchNumber = 0;
        boolean init = true;
        ProgressLogger pgEpoch = new ProgressLogger(LOG);
        pgEpoch.itemsName = "epoch";
        pgEpoch.expectedUpdates = numEpochs;
        pgEpoch.start();
        double bestScore = Double.MAX_VALUE;
        LocalFileModelSaver saver = new LocalFileModelSaver(directory);
        int iter = 0;
        Map<Integer, Double> scoreMap=new HashMap<Integer,Double>();
        for (int epoch = 0; epoch < numEpochs; epoch++) {
            ProgressLogger pg = new ProgressLogger(LOG);
            pg.itemsName = "mini-batch";

            pg.expectedUpdates = async.totalExamples() / miniBatchSize; // one iteration processes miniBatchIterator elements.
            pg.start();
            int lastIter = 0;
            while (async.hasNext()) {
                // trigger bugs early:
//                printSample(miniBatchSize, exampleLength, nSamplesToGenerate, nCharactersToSample, generationInitialization, rng, directory, iter, net, miniBatchNumber);

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
                    System.out.println(ds.toString());
                    System.err.println("Aborting because NaN was generated for score.");
                    System.exit(1);
                }
                iter++;
                if (score < bestScore * 0.95 && iter > (lastIter + MIN_ITERATION_BETWEEN_BEST_MODEL)) {
                    bestScore = score;
                    saver.saveBestModel(net, score);
                    System.out.println("Saving best score model.. score=" + bestScore);
                    lastIter = iter;
                }
                scoreMap.put(iter,bestScore);
            }
            //save latest after the end of an epoch:
            saver.saveLatestModel(net, net.score());
            saveModel(saver, directory, Integer.toString(epoch) + "-", net);

            pg.stop();
            pgEpoch.update();
            async.reset();    //Reset iterator for another epoch
        }
        pgEpoch.stop();

        return new EarlyStoppingResult<MultiLayerNetwork>(EarlyStoppingResult.TerminationReason.EpochTerminationCondition,
         "not early stopping",scoreMap, numEpochs,bestScore, numEpochs, net);
    }
}
