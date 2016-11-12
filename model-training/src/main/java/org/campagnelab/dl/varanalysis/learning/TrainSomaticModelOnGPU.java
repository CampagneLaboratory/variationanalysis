package org.campagnelab.dl.varanalysis.learning;

import it.unimi.dsi.logging.ProgressLogger;
import org.campagnelab.dl.varanalysis.learning.iterators.NamedCachingDataSetIterator;
import org.campagnelab.dl.varanalysis.learning.iterators.NamedDataSetIterator;
import org.campagnelab.dl.varanalysis.learning.models.ModelSaver;
import org.campagnelab.dl.varanalysis.learning.performance.MeasurePerformance;
import org.deeplearning4j.earlystopping.EarlyStoppingResult;
import org.deeplearning4j.nn.conf.MultiLayerConfiguration;
import org.deeplearning4j.nn.multilayer.MultiLayerNetwork;
import org.deeplearning4j.parallelism.ParallelWrapper;
import org.nd4j.linalg.api.buffer.DataBuffer;
import org.nd4j.linalg.api.buffer.util.DataTypeUtil;
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

//import org.nd4j.jita.conf.CudaEnvironment;

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

    private String validationDatasetFilename = null;


    @Override
    public void execute() {
        if (args().trainingSets.size() == 0) {
            System.out.println("Please add at least one training set to the args().");
            return;
        }
  /*  CudaEnvironment.getInstance().getConfiguration()
                .enableDebug(false)
                .allowMultiGPU(true)
                .setMaximumGridSize(512)
                .setMaximumBlockSize(512)
                .setMaximumDeviceCacheableLength(1024 * 1024 * 1024L)
                .setMaximumDeviceCache(8L * 1024 * 1024 * 1024L)
                .setMaximumHostCacheableLength(1024 * 1024 * 1024L)
                .setMaximumHostCache(8L * 1024 * 1024 * 1024L)
                // cross-device access is used for faster model averaging over pcie
                .allowCrossDeviceAccess(true);*/
        if ("FP16".equals(args().precision)) {
            precision = ParameterPrecision.FP16;
            System.out.println("Parameter precision set to FP16.");
        }
        super.execute();
    }

    public static void main(String[] args) throws IOException {

        TrainSomaticModelOnGPU tool = new TrainSomaticModelOnGPU();
        SomaticTrainingArguments arguments = tool.createArguments();
        tool.parseArguments(args, "TrainSomaticModelOnGPU", arguments);
        if ("FP16".equals(tool.args().precision)) {
            DataTypeUtil.setDTypeForContext(DataBuffer.Type.HALF);
        }
        tool.execute();
        tool.writeModelingConditions(arguments);
        System.err.println("Allow Multi-GPU");
    }


    @Override
    protected NamedDataSetIterator decorateIterator(NamedDataSetIterator iterator) {
        return new NamedCachingDataSetIterator(new CachingDataSetIterator(iterator, new InMemoryDataSetCache()), "<no basename defined>");
    }

    @Override
    protected EarlyStoppingResult<MultiLayerNetwork> train(MultiLayerConfiguration conf, DataSetIterator async) throws IOException {
        validationDatasetFilename = args().validationSet;
        //check validation file for error
        if (!(new File(validationDatasetFilename).exists())) {
            throw new IOException("Validation file not found! " + validationDatasetFilename);
        }
        perf = new MeasurePerformance(args().numValidation, validationDatasetFilename, args().miniBatchSize, featureCalculator, labelMapper);

        ParallelWrapper wrapper = new ParallelWrapper.Builder(net)
                .prefetchBuffer(args().miniBatchSize)
                .workers(4)
                .averagingFrequency(1)
                .reportScoreAfterAveraging(false)
                .useLegacyAveraging(false)
                .build();

        //Do training, and then generate and print samples from network
        int miniBatchNumber = 0;
        boolean init = true;
        ProgressLogger pgEpoch = new ProgressLogger(LOG);
        pgEpoch.displayLocalSpeed = true;
        pgEpoch.itemsName = "epoch";
        pgEpoch.expectedUpdates = args().maxEpochs;
        pgEpoch.start();
        bestScore = Double.MAX_VALUE;
        ModelSaver saver = new ModelSaver(directory);
        int numExamplesUsed = 0;

        Map<Integer, Double> scoreMap = new HashMap<Integer, Double>();
        double bestAUC = 0;
        int notImproved = 0;
        int iter = 0;
        int epoch;
        assert async.resetSupported(): "Iterator must support reset.";
        for (epoch = 0; epoch < args().maxEpochs; epoch++) {

            wrapper.fit(async);
            pgEpoch.update();

            double score = net.score();
            scoreMap.put(epoch, score);
            bestScore = Math.min(score, bestScore);

            writeBestScoreFile();
            async.reset();

            saver.saveLatestModel(net, net.score());

            writeProperties(this);
            writeBestScoreFile();
            if (epoch % args().validateEvery == 0) {
                double auc = estimateTestSetPerf(epoch, iter);
                performanceLogger.log("epochs", numExamplesUsed, epoch, score, auc);
                if (auc > bestAUC) {
                    saver.saveModel(net, "bestAUC", auc);
                    bestAUC = auc;
                    writeBestAUC(bestAUC);
                    performanceLogger.log("bestAUC", numExamplesUsed, epoch, bestScore, bestAUC);
                    notImproved = 0;
                } else {
                    notImproved++;
                }
                if (notImproved > args().stopWhenEpochsWithoutImprovement) {
                    // we have not improved after earlyStopCondition epoch, time to stop.
                    break;
                }
                System.out.printf("epoch %d auc=%g%n", epoch, auc);
            }

            numExamplesUsed += async.totalExamples();
            performanceLogger.write();
        }

        pgEpoch.stop();
        //TODO enable with 0.6.0+
        //     wrapper.shutdown();
        return new EarlyStoppingResult<MultiLayerNetwork>(EarlyStoppingResult.TerminationReason.EpochTerminationCondition,
                "not early stopping", scoreMap, performanceLogger.getBestEpoch("bestAUC"), bestScore, args().maxEpochs, net);
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
        double auc = perf.estimateAUC(net);

        System.out.printf("Epoch %d Iteration %d AUC=%f %n", epoch, iter, auc);

        return auc;
    }

    @Override
    public SomaticTrainingArguments createArguments() {
        return new SomaticTrainingArguments();
    }
}
