package org.campagnelab.dl.varanalysis.learning;

import it.unimi.dsi.logging.ProgressLogger;
import org.campagnelab.dl.model.utils.mappers.FeatureMapperV15;
import org.campagnelab.dl.varanalysis.learning.parallel.ParallelWrapper;
import org.deeplearning4j.earlystopping.EarlyStoppingResult;
import org.deeplearning4j.earlystopping.saver.LocalFileModelSaver;
import org.deeplearning4j.nn.conf.MultiLayerConfiguration;
import org.deeplearning4j.nn.multilayer.MultiLayerNetwork;

import org.nd4j.linalg.api.ndarray.INDArray;
import org.nd4j.linalg.dataset.DataSet;
import org.nd4j.linalg.dataset.api.iterator.DataSetIterator;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

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


    public static void main(String[] args) throws IOException {
        // uncomment the following line when running on a machine with multiple GPUs:
        //  org.nd4j.jita.conf.CudaEnvironment.getInstance().getConfiguration().allowMultiGPU(true);
        TrainingArguments arguments = parseArguments(args, "TrainSomaticModelOnGPU");
        TrainSomaticModelOnGPU trainer = new TrainSomaticModelOnGPU(arguments);
        trainer.execute(new FeatureMapperV15(), arguments.getTrainingSets(), arguments.miniBatchSize);
    }


    @Override
    protected EarlyStoppingResult<MultiLayerNetwork> train(MultiLayerConfiguration conf, DataSetIterator async) throws IOException {

        ParallelWrapper wrapper = new ParallelWrapper.Builder(net)
                .prefetchBuffer(arguments.miniBatchSize)
                .workers(arguments.miniBatchSize)
                .averagingFrequency(3)
                .build();

        //Do training, and then generate and print samples from network
        int miniBatchNumber = 0;
        boolean init = true;
        ProgressLogger pgEpoch = new ProgressLogger(LOG);
        pgEpoch.itemsName = "epoch";
        pgEpoch.expectedUpdates = arguments.maxEpochs;
        pgEpoch.start();
        bestScore = Double.MAX_VALUE;
        LocalFileModelSaver saver = new LocalFileModelSaver(directory);

        Map<Integer, Double> scoreMap = new HashMap<Integer, Double>();
        for (int epoch = 0; epoch < arguments.maxEpochs; epoch++) {

            wrapper.fit(async);
            pgEpoch.update();
            async.reset();
            double score = net.score();
            scoreMap.put(epoch, score);
            bestScore = Math.min(score, bestScore);
            saveModel(saver, directory, String.format("%d-", epoch), net);

            writeBestScoreFile();
        }

        pgEpoch.stop();

        return new EarlyStoppingResult<MultiLayerNetwork>(EarlyStoppingResult.TerminationReason.EpochTerminationCondition,
                "not early stopping", scoreMap, arguments.maxEpochs, bestScore, arguments.maxEpochs, net);
    }
}
