package org.campagnelab.dl.varanalysis.learning;

import it.unimi.dsi.logging.ProgressLogger;
import org.campagnelab.dl.model.utils.mappers.FeatureMapperV18;
import org.campagnelab.dl.model.utils.mappers.trio.FeatureMapperV18Trio;
import org.campagnelab.dl.varanalysis.learning.models.ModelSaver;
import org.deeplearning4j.earlystopping.EarlyStoppingResult;
import org.deeplearning4j.nn.conf.MultiLayerConfiguration;
import org.deeplearning4j.nn.multilayer.MultiLayerNetwork;
import org.deeplearning4j.parallelism.ParallelWrapper;
import org.nd4j.linalg.dataset.api.iterator.CachingDataSetIterator;
import org.nd4j.linalg.dataset.api.iterator.DataSetIterator;
import org.nd4j.linalg.dataset.api.iterator.cache.InMemoryDataSetCache;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.HashMap;
import java.util.Map;

/**
 * Train a neural network to predict mutations.
 * <p>
 * Created by fac2003 on 5/21/16.
 *
 * @author Fabien Campagne
 */
public class TrainSomaticModelOnGPU extends SomaticTrainer {
    public static final int MIN_ITERATION_BETWEEN_BEST_MODEL = 1000;
    static private Logger LOG = LoggerFactory.getLogger(TrainSomaticModelOnGPU.class);

    public TrainSomaticModelOnGPU(TrainingArguments arguments) {
        super(arguments);
    }

    private String validationDatasetFilename = null;


    public static void main(String[] args) throws IOException {

        // uncomment the following line when running on a machine with multiple GPUs:
        //  org.nd4j.jita.conf.CudaEnvironment.getInstance().getConfiguration().allowMultiGPU(true);
        TrainingArguments arguments = parseArguments(args, "TrainSomaticModelOnGPU");
        TrainSomaticModelOnGPU trainer = new TrainSomaticModelOnGPU(arguments);
        if (arguments.isTrio) {
            trainer.execute(new FeatureMapperV18Trio(), arguments.getTrainingSets(), arguments.miniBatchSize);
        } else {
            trainer.execute(new FeatureMapperV18(), arguments.getTrainingSets(), arguments.miniBatchSize);
        }
    }

    @Override
    protected DataSetIterator decorateIterator(DataSetIterator iterator) {
        return new CachingDataSetIterator(iterator, new InMemoryDataSetCache());
    }

    @Override
    protected EarlyStoppingResult<MultiLayerNetwork> train(MultiLayerConfiguration conf, DataSetIterator async) throws IOException {
        validationDatasetFilename = arguments.validationSet;
        //check validation file for error
        if (!(new File(validationDatasetFilename).exists())) {
            throw new IOException("Validation file not found! " + validationDatasetFilename);
        }
        perf = new MeasurePerformance(arguments.numValidation, validationDatasetFilename);

        ParallelWrapper wrapper = new ParallelWrapper.Builder(net)
                .prefetchBuffer(arguments.miniBatchSize)
                .workers(4)
                .averagingFrequency(1)
                .reportScoreAfterAveraging(false)
                .useLegacyAveraging(false)
                .build();

        //Do training, and then generate and print samples from network
        int miniBatchNumber = 0;
        boolean init = true;
        ProgressLogger pgEpoch = new ProgressLogger(LOG);
        pgEpoch.itemsName = "epoch";
        pgEpoch.expectedUpdates = arguments.maxEpochs;
        pgEpoch.start();
        bestScore = Double.MAX_VALUE;
        ModelSaver saver = new ModelSaver(directory);
        int numExamplesUsed = 0;

        Map<Integer, Double> scoreMap = new HashMap<Integer, Double>();
        double bestAUC = 0;
        int notImproved = 0;
        int iter = 0;
        for (int epoch = 0; epoch < arguments.maxEpochs; epoch++) {

            wrapper.fit(async);
            pgEpoch.update();

            double score = net.score();
            scoreMap.put(epoch, score);
            bestScore = Math.min(score, bestScore);

            async.reset();
            writeBestScoreFile();
            if (async.resetSupported()) {
                async.reset();
            }

            writeProperties(this);
            writeBestScoreFile();
            double auc = estimateTestSetPerf(epoch, iter);
            performanceLogger.log("epochs", numExamplesUsed, epoch, Double.NaN, auc);
            if (auc > bestAUC) {
                saver.saveModel(net, "bestAUC", auc);
                bestAUC = auc;
                writeBestAUC(bestAUC);
                performanceLogger.log("bestAUC", numExamplesUsed, epoch, bestScore, bestAUC);
                notImproved = 0;
            } else {
                notImproved++;
            }
            if (notImproved > arguments.stopWhenEpochsWithoutImprovement) {
                // we have not improved after earlyStopCondition epoch, time to stop.
                break;
            }
            System.out.printf("epoch %d auc=%g%n", epoch, auc);
            numExamplesUsed += arguments.numTraining;
        }

        pgEpoch.stop();

        return new EarlyStoppingResult<MultiLayerNetwork>(EarlyStoppingResult.TerminationReason.EpochTerminationCondition,
                "not early stopping", scoreMap, arguments.maxEpochs, bestScore, arguments.maxEpochs, net);
    }

    private void writeBestAUC(double bestAUC) {
        try {
            FileWriter scoreWriter = new FileWriter(directory + "/bestAUC");
            scoreWriter.append(Double.toString(bestAUC));
            scoreWriter.close();
        } catch (IOException e) {

        }

    }

    MeasurePerformance perf;

    protected double estimateTestSetPerf(int epoch, int iter) throws IOException {
        if (validationDatasetFilename == null) return 0;
        double auc = perf.estimateAUC(featureCalculator, net);

        System.out.printf("Epoch %d Iteration %d AUC=%f %n", epoch, iter, auc);

        return auc;
    }
}
