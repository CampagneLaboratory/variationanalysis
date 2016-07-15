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
import org.campagnelab.dl.varanalysis.util.HitBoundedPriorityQueue;
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
    public static final int MAX_ERRORS_KEPT = 1;
    static private Logger LOG = LoggerFactory.getLogger(TrainSomaticModel.class);
    private String validationDatasetFilename = null;


    public static void main(String[] args) throws IOException {

        TrainSomaticModel trainer = new TrainSomaticModel();
        if (args.length < 1) {
            System.err.println("usage: DetectMutations <input-training-file+>");
        }
        trainer.execute(new FeatureMapperV14(), args, 32);

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
        bestScore = Double.MAX_VALUE;
        LocalFileModelSaver saver = new LocalFileModelSaver(directory);
        int iter = 0;
        Map<Integer, Double> scoreMap = new HashMap<Integer, Double>();

        for (int epoch = 0; epoch < numEpochs; epoch++) {
            ProgressLogger pg = new ProgressLogger(LOG);
            pg.itemsName = "mini-batch";

            pg.expectedUpdates = async.totalExamples() / miniBatchSize; // one iteration processes miniBatchIterator elements.
            pg.start();
            int lastIter = 0;
            while (async.hasNext()) {

                DataSet ds = async.next();
                if (numLabels(ds.getLabels()) != 2) {
                    System.out.println("There should be two labels in the miniBatch");
                }

                DataSet enriched=new DataSet();
                // fit the net:
                net.fit(ds);

                INDArray predictedLabels = net.output(ds.getFeatures(), false);
                keepWorseErrors(ds, predictedLabels, ds.getLabels());


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
                    estimateTestSetPerf(epoch, iter);
                }
                scoreMap.put(iter, bestScore);
                queue.clear();
            }
            //save latest after the end of an epoch:
            saver.saveLatestModel(net, net.score());
            saveModel(saver, directory, Integer.toString(epoch) + "-", net);
            writeProperties(featureCalculator);
            writeBestScoreFile();
            estimateTestSetPerf(epoch, iter);
            pg.stop();
            pgEpoch.update();
            async.reset();    //Reset iterator for another epoch
        }
        pgEpoch.stop();

        return new EarlyStoppingResult<MultiLayerNetwork>(EarlyStoppingResult.TerminationReason.EpochTerminationCondition,
                "not early stopping", scoreMap, numEpochs, bestScore, numEpochs, net);
    }

    HitBoundedPriorityQueue queue = new HitBoundedPriorityQueue(MAX_ERRORS_KEPT);

    private void keepWorseErrors( DataSet minibatch, INDArray predictedLabels, INDArray labels) {
        int size=minibatch.numExamples();
        for (int exampleIndex = 0; exampleIndex < size; exampleIndex++) {
            if (isWrongPrediction(exampleIndex, predictedLabels, labels)) {
                float wrongness = calculateWrongness(exampleIndex, predictedLabels, labels);
                queue.enqueue(wrongness, minibatch.getFeatures(), labels.getRow(exampleIndex));
            //    System.out.println("largest error so far: "+ queue.first());
            }
        }
    }

    private float calculateWrongness(int exampleIndex, INDArray predictedLabels, INDArray labels) {
        final int positiveLabelMutated = 0;
        final int negativeLabel = 1;
        final boolean predictedPositive = predictedLabels.getDouble(exampleIndex, positiveLabelMutated) > predictedLabels.getDouble(exampleIndex, negativeLabel);
        if (predictedPositive) {
            return predictedLabels.getFloat(exampleIndex, positiveLabelMutated);
        } else{
            return predictedLabels.getFloat(exampleIndex, negativeLabel);
        }

    }

    private boolean isWrongPrediction(int exampleIndex, INDArray predictedLabels, INDArray labels) {
        final int positiveLabelMutated = 0;
        final int negativeLabel = 1;
        final boolean predictedPositive = predictedLabels.getDouble(exampleIndex, positiveLabelMutated) > predictedLabels.getDouble(exampleIndex, negativeLabel);
        final double trueLabelPositive = labels.getDouble(exampleIndex, positiveLabelMutated);
        return predictedPositive && trueLabelPositive ==0 ||
                !predictedPositive && trueLabelPositive >0;
    }

    private void estimateTestSetPerf(int epoch, int iter) throws IOException {
        if (validationDatasetFilename == null) return;
        MeasurePerformance perf = new MeasurePerformance(10000);
        validationDatasetFilename = "/data/mutated-MHFC-13-CTL_B_NK.parquet";
        double auc = perf.estimateAUC(featureCalculator, net, validationDatasetFilename);
        System.out.printf("Epoch %d Iteration %d AUC=%f%n", epoch, iter, auc);
    }
}
