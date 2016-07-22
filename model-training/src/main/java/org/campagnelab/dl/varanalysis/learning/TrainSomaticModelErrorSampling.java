package org.campagnelab.dl.varanalysis.learning;

import it.unimi.dsi.logging.ProgressLogger;
import org.campagnelab.dl.model.utils.mappers.FeatureMapperV16;
import org.campagnelab.dl.varanalysis.learning.io.ModelSaver;
import org.campagnelab.dl.varanalysis.learning.iterators.BaseInformationConcatIterator;
import org.campagnelab.dl.varanalysis.learning.iterators.SamplingIterator;
import org.campagnelab.dl.varanalysis.util.ErrorRecord;
import org.deeplearning4j.earlystopping.EarlyStoppingResult;
import org.deeplearning4j.nn.conf.MultiLayerConfiguration;
import org.deeplearning4j.nn.multilayer.MultiLayerNetwork;
import org.nd4j.linalg.api.ndarray.INDArray;
import org.nd4j.linalg.dataset.DataSet;
import org.nd4j.linalg.dataset.api.iterator.DataSetIterator;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.io.FileWriter;
import java.io.IOException;
import java.util.HashMap;
import java.util.Map;

/**
 * Train a neural network to predict mutations. Use error sampling to focus learning on parts of the training
 * set where errors are still being made. This should result in lowering the probability assigned to frequent
 * types of errors.
 * <p>
 *
 * @author Fabien Campagne
 */
public class TrainSomaticModelErrorSampling extends SomaticTrainer {

    private static final int MAX_EPOCHS_WITHOUT_IMPROVEMENT = 10;
    static private Logger LOG = LoggerFactory.getLogger(TrainSomaticModelErrorSampling.class);
    private String validationDatasetFilename = null;

    /**
     * Error enrichment support.
     **/

    private final boolean ERROR_SAMPLING = true;


    public static void main(String[] args) throws IOException {

        TrainSomaticModelErrorSampling trainer = new TrainSomaticModelErrorSampling();
        if (args.length < 1) {
            System.err.println("usage: DetectMutations <input-training-file+>");
        }
        trainer.execute(new FeatureMapperV16(), args, 32);

    }

    SamplingIterator samplingIterator = null;

    @Override
    protected DataSetIterator decorateIterator(BaseInformationConcatIterator iterator) {
        samplingIterator = new SamplingIterator(iterator, seed);
        return samplingIterator;
    }

    @Override
    protected EarlyStoppingResult<MultiLayerNetwork> train(MultiLayerConfiguration conf, DataSetIterator async) throws IOException {
        //Do training, and then generate and print samples from network
        int miniBatchNumber = 0;
        boolean init = true;
        ProgressLogger pgEpoch = new ProgressLogger(LOG);
        pgEpoch.itemsName = "epoch";
        pgEpoch.expectedUpdates = numEpochs;
        pgEpoch.displayLocalSpeed = true;
        pgEpoch.start();
        bestScore = Double.MAX_VALUE;
        ModelSaver saver = new ModelSaver(directory);
        int iter = 0;
        int numEpochSinceImprovement = 0;
        Map<Integer, Double> scoreMap = new HashMap<Integer, Double>();
        System.out.println("ERROR_SAMPLING=" + ERROR_SAMPLING);
        double bestAUC = 0.5;
        double finalAUC = 0.5;
        for (int epoch = 0; epoch < numEpochs; epoch++) {
            ProgressLogger pg = new ProgressLogger(LOG);
            pg.itemsName = "mini-batch";
            pg.displayLocalSpeed = true;
            pg.expectedUpdates = async.totalExamples() / miniBatchSize; // one iteration processes miniBatchIterator elements.
            pg.start();
            int maxProcessPerEpoch = Integer.MAX_VALUE;
            while (async.hasNext()) {

                DataSet ds = async.next();
                if (numLabels(ds.getLabels()) != 2) {
                    // just be chance, some minibatches will have zero or only one label.
                    // Since we cannot learn from them, we skip them.
                    continue;
                }

                // fit the net:
                net.fit(ds);

                INDArray predictedLabels = net.output(ds.getFeatures(), false);
                updateProbabilities(predictedLabels, ds.getLabels());

                pg.update();
                iter++;

                if (iter > maxProcessPerEpoch) break;

            }

            pg.stop();
            pgEpoch.update();
            async.reset();    //Reset iterator for another epoch

            double auc = estimateTestSetPerf(epoch, iter);
            if (auc > bestAUC) {
                saver.saveModel(net, "bestAUC", auc);
                bestAUC = auc;
                writeBestAUC(bestAUC);
                writeProperties(this);
                numEpochSinceImprovement = 0;
            } else {
                numEpochSinceImprovement++;
                if (numEpochSinceImprovement > MAX_EPOCHS_WITHOUT_IMPROVEMENT) {
                    System.out.printf("AUC did not improve after %d epoch. Early stop.",
                            MAX_EPOCHS_WITHOUT_IMPROVEMENT);
                    finalAUC = auc;
                    break;
                }
            }

        }
        pgEpoch.stop();
        saver.saveModel(net, "final", finalAUC);
        return new EarlyStoppingResult<MultiLayerNetwork>(EarlyStoppingResult.TerminationReason.EpochTerminationCondition,
                "not early stopping", scoreMap, numEpochs, bestScore, numEpochs, net);
    }

    /**
     * Update the iterator with the probability that each record of the minibatch should be reused.
     *
     * @param predictedLabels
     * @param labels
     */
    private void updateProbabilities(INDArray predictedLabels, INDArray labels) {
        for (int exampleIndex = 0; exampleIndex < predictedLabels.size(0); exampleIndex++) {

            final float pOfWrongLabel = ErrorRecord.calculateWrongness(exampleIndex, predictedLabels, labels);
            final boolean wrongPrediction = ErrorRecord.isWrongPrediction(exampleIndex, predictedLabels, labels);
            float p = wrongPrediction ? pOfWrongLabel :
                    1 - pOfWrongLabel;
            /*if (!wrongPrediction) {
                p=0.05f;
            }*/
            final float previousP = samplingIterator.getProbability(exampleIndex);
            samplingIterator.setSamplingProbability(wrongPrediction, exampleIndex, p * previousP);

        }
    }


    private void writeBestAUC(double bestAUC) {
        try {
            FileWriter scoreWriter = new FileWriter(directory + "/bestAUC");
            scoreWriter.append(Double.toString(bestAUC));
            scoreWriter.close();
        } catch (IOException e) {

        }

    }

    @Override
    public void appendProperties(ModelPropertiesHelper helper) {
        super.appendProperties(helper);
        helper.put("ErrorEnrichment.active", ERROR_SAMPLING);
    }


    private double estimateTestSetPerf(int epoch, int iter) throws IOException {
        validationDatasetFilename = "/data/no-threshold/validation-2/m-r-MHFC-63-CTL_B_NK_VN.parquet";
        if (validationDatasetFilename == null) return 0;
        MeasurePerformance perf = new MeasurePerformance(10000);
        double auc = perf.estimateAUC(featureCalculator, net, validationDatasetFilename);
        //   System.out.printf("Average sampling P=%f%n", samplingIterator.getAverageSamplingP());

        System.out.printf("Epoch %d Iteration %d AUC=%f Average sampling P=%f %%skipped=%f %n ", epoch, iter, auc, samplingIterator.getAverageSamplingP(),
                samplingIterator.percentSkipped());

        return auc;
    }
}
