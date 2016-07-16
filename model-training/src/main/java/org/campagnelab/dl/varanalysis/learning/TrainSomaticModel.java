package org.campagnelab.dl.varanalysis.learning;

import it.unimi.dsi.fastutil.objects.ObjectArrayList;
import it.unimi.dsi.logging.ProgressLogger;
import org.campagnelab.dl.model.utils.mappers.*;
import org.campagnelab.dl.varanalysis.util.ErrorRecord;
import org.campagnelab.dl.varanalysis.util.HitBoundedPriorityQueue;
import org.deeplearning4j.earlystopping.EarlyStoppingResult;
import org.deeplearning4j.earlystopping.saver.LocalFileModelSaver;
import org.deeplearning4j.nn.conf.MultiLayerConfiguration;
import org.deeplearning4j.nn.multilayer.MultiLayerNetwork;
import org.nd4j.linalg.api.ndarray.INDArray;
import org.nd4j.linalg.dataset.DataSet;
import org.nd4j.linalg.dataset.api.iterator.DataSetIterator;
import org.nd4j.linalg.factory.Nd4j;
import org.nd4j.linalg.indexing.PointIndex;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.io.*;
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
    private String validationDatasetFilename = null;

    /**
     * Error enrichment support.
     **/
    public static final int MAX_ERRORS_KEPT = 25;
    private final boolean ERROR_ENRICHMENT = true;
    private final int NUM_ERRORS_ADDED = 3;
    private HitBoundedPriorityQueue queue = new HitBoundedPriorityQueue(MAX_ERRORS_KEPT);

    public static void main(String[] args) throws IOException {

        TrainSomaticModel trainer = new TrainSomaticModel();
        if (args.length < 1) {
            System.err.println("usage: DetectMutations <input-training-file+>");
        }
        trainer.execute(new FeatureMapperV13(), args, 32);

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
        System.out.println("ERROR_ENRICHMENT=" + ERROR_ENRICHMENT);
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

                ds = enrichWithErrors(ds);
                // fit the net:
                net.fit(ds);

                INDArray predictedLabels = net.output(ds.getFeatures(), false);
                keepWorseErrors(ds, predictedLabels, ds.getLabels());
                pg.update();
                queue.updateWrongness(net);
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

            }
            //save latest after the end of an epoch:
            saver.saveLatestModel(net, net.score());
            saveModel(saver, directory, Integer.toString(epoch) + "-", net);
            writeProperties(this);
            writeBestScoreFile();
            estimateTestSetPerf(epoch, iter);
            pg.stop();
            pgEpoch.update();
            queue.clear();
            async.reset();    //Reset iterator for another epoch
        }
        pgEpoch.stop();

        return new EarlyStoppingResult<MultiLayerNetwork>(EarlyStoppingResult.TerminationReason.EpochTerminationCondition,
                "not early stopping", scoreMap, numEpochs, bestScore, numEpochs, net);
    }

    @Override
    public void appendProperties(ModelPropertiesHelper helper) {
        super.appendProperties(helper);
    helper.put("ErrorEnrichment.active",ERROR_ENRICHMENT);
    helper.put("ErrorEnrichment.MAX_ERRORS_KEPT",MAX_ERRORS_KEPT);
    helper.put("ErrorEnrichment.NUM_ERRORS_ADDED",NUM_ERRORS_ADDED);
    }

    DataSet[] array = new DataSet[2];

    private DataSet enrichWithErrors(DataSet ds) {
        if (!ERROR_ENRICHMENT) {
            return ds;
        }
        if (queue.isEmpty()) {
            // no errors were collected yet. Return the un-enriched dataset.
            return ds;
        }

        int size = this.NUM_ERRORS_ADDED;
        INDArray inputs = Nd4j.zeros(size, featureCalculator.numberOfFeatures());
        INDArray labels = Nd4j.zeros(size, labelMapper.numberOfLabels());
        int i = 0;

        for (ErrorRecord errorRecord : queue.getRandomSample(size)) {
            // we are going to call nextRecord directly, without checking hasNextRecord, because we have
            // determined how many times we can call (in size). We should get the exception if we were
            // wrong in our estimate of size.

            // fill in features and labels for a given record i:
            Nd4j.copy(errorRecord.features.get(new PointIndex(0)), inputs.get(new PointIndex(i)));
            Nd4j.copy(errorRecord.label, labels.get(new PointIndex(i)));
            i++;

        }
        DataSet errorDataSet = new DataSet(inputs, labels);
        array[0] = ds;
        array[1] = errorDataSet;
        final DataSet enrichedDataset = DataSet.merge(ObjectArrayList.wrap(array));
        return enrichedDataset;
    }


    private void keepWorseErrors(DataSet minibatch, INDArray predictedLabels, INDArray labels) {
        if (!ERROR_ENRICHMENT) {
            return;
        }
        int size = minibatch.numExamples();
        for (int exampleIndex = 0; exampleIndex < size; exampleIndex++) {
            if (isWrongPrediction(exampleIndex, predictedLabels, labels)) {
                float wrongness = calculateWrongness(exampleIndex, predictedLabels, labels);
                queue.enqueue(wrongness, minibatch.getFeatures(), labels.getRow(exampleIndex));
                //    System.out.println("largest error so far: "+ queue.first());
            }
        }
    }

    private float calculateWrongness(int exampleIndex, INDArray predictedLabels, INDArray labels) {

        return ErrorRecord.calculateWrongness(exampleIndex, predictedLabels, labels);

    }

    private boolean isWrongPrediction(int exampleIndex, INDArray predictedLabels, INDArray labels) {
        if (!ERROR_ENRICHMENT) {
            return false;
        }
        return ErrorRecord.isWrongPrediction(exampleIndex, predictedLabels, labels);

    }

    private void estimateTestSetPerf(int epoch, int iter) throws IOException {
        validationDatasetFilename = "/data/NN/validation/mutated-MHFC-53-CTL_B_NK.parquet";
        if (validationDatasetFilename == null) return;
        MeasurePerformance perf = new MeasurePerformance(10000);
        double auc = perf.estimateAUC(featureCalculator, net, validationDatasetFilename);
        float minWrongness=queue.getMinWrongness();
        float maxWrongness=queue.getMaxWrongness();
        float meanWrongness=queue.getMeanWrongness();
        System.out.printf("Epoch %d Iteration %d AUC=%f wrongness: %f-%f mean: %f #examples: %d %n", epoch, iter, auc,
                minWrongness, maxWrongness,meanWrongness,
                queue.size());
    }
}
