package org.campagnelab.dl.varanalysis.learning.calibrate;

import it.unimi.dsi.fastutil.floats.FloatArrayList;
import it.unimi.dsi.logging.ProgressLogger;
import org.campagnelab.dl.model.utils.ProtoPredictor;
import org.campagnelab.dl.model.utils.mappers.FeatureMapper;
import org.campagnelab.dl.model.utils.models.ModelLoader;
import org.campagnelab.dl.varanalysis.learning.architecture.CalibrationAssembler;
import org.campagnelab.dl.varanalysis.learning.architecture.NeuralNetAssembler;
import org.campagnelab.dl.varanalysis.learning.models.ModelSaver;
import org.campagnelab.dl.varanalysis.protobuf.BaseInformationRecords;
import org.campagnelab.dl.varanalysis.storage.RecordReader;
import org.deeplearning4j.nn.api.Layer;
import org.deeplearning4j.nn.conf.LearningRatePolicy;
import org.deeplearning4j.nn.conf.MultiLayerConfiguration;
import org.deeplearning4j.nn.multilayer.MultiLayerNetwork;
import org.deeplearning4j.nn.weights.WeightInit;
import org.deeplearning4j.optimize.listeners.ScoreIterationListener;
import org.nd4j.linalg.api.ndarray.INDArray;
import org.nd4j.linalg.dataset.DataSet;
import org.nd4j.linalg.factory.Nd4j;
import org.nd4j.linalg.lossfunctions.LossFunctions;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.io.IOException;

/**
 * Created by fac2003 on 7/17/16.
 */
public class CalibrateModelOutput {
    private static final int MIN_ITERATION_BETWEEN_BEST_MODEL = 10000;
    static private Logger LOG = LoggerFactory.getLogger(CalibrateModelOutput.class);
    private int seed = 232323;
    private double learningRate = 0.1;
    private int miniBatchSize = 32;
    private int numEpoch = 3;


    public static void main(String[] args) throws IOException {
        if (args.length < 2) {
            System.err.println("Usage: modelDirectory modelPrefix calibrationDatafile");
            System.exit(1);
        }
        String calibrationDatafile = args[2];
        String modelDirectory = args[0];
        String modelPrefix = args[1];
        CalibrateModelOutput calibrator = new CalibrateModelOutput(modelDirectory, modelPrefix, calibrationDatafile);
        calibrator.calibrate();
    }

    private void calibrate() throws IOException {
        ModelLoader loader = new ModelLoader(modelDirectory);

        MultiLayerNetwork model = null;
        try {
            model = loader.loadModel(modelPrefix);
        } catch (IOException e) {
            System.err.printf("Unable to load model with modelPath=%s and prefix=%s%n", modelDirectory, modelPrefix);
            System.exit(1);
        }

        RecordReader reader = null;
        try {
            reader = new RecordReader(calibrationDatafile);
        } catch (IOException e) {
            System.err.printf("Unable to load calibration file %s%n", calibrationDatafile);
            System.exit(1);
        }


        //DataSet ds = baseIter.next();
//set up logger

        MultiLayerNetwork calibrationModel;
        FeatureMapper modelFeatureMapper = loader.loadFeatureMapper();

        int featureNumber = getModelActivationNumber(model, modelFeatureMapper);
        int numInputs = featureNumber;
        int numOutputs = 1;
        int numHiddenNodes = numInputs;
        NeuralNetAssembler assembler = new CalibrationAssembler();
        assembler.setSeed(seed);
        assembler.setLearningRate(learningRate);
        assembler.setNumHiddenNodes(numHiddenNodes);
        assembler.setNumInputs(numInputs);
        assembler.setNumOutputs(numOutputs);
        assembler.setLossFunction(LossFunctions.LossFunction.RMSE_XENT); // RMSE for regression.
        assembler.setRegularization(true);
        assembler.setRegularizationRate(1e-6);

        assembler.setWeightInitialization(WeightInit.RELU);
        assembler.setLearningRatePolicy(LearningRatePolicy.Score);
        MultiLayerConfiguration conf = assembler.createNetwork();
        calibrationModel = new MultiLayerNetwork(conf);
        calibrationModel.init();
        calibrationModel.setListeners(new ScoreIterationListener(100));

        int iter = 0;
        int lastIter = 0;
        DataSet miniBatch = null;
        int currentMiniBatchSize = 0;
        ModelSaver saver = new ModelSaver(modelDirectory);
        double score = Float.MAX_VALUE;
        double bestScore = Float.MAX_VALUE - 10;
        for (int epoch = 0; epoch < numEpoch; epoch++) {
            ProgressLogger pgReadWrite = new ProgressLogger(LOG);
            pgReadWrite.itemsName = "example";
            pgReadWrite.expectedUpdates = reader.getTotalRecords();

            pgReadWrite.start();
            for (BaseInformationRecords.BaseInformation record : reader) {

                if (currentMiniBatchSize == 0) {

                    // create a new empty mini-batch. Done at the beginning and every time we are done with a
                    // completed mini-batch:
                    INDArray inputs = Nd4j.zeros(miniBatchSize, featureNumber);
                    INDArray labels = Nd4j.zeros(miniBatchSize, 1);
                    miniBatch = new DataSet(inputs, labels);
                }
                currentMiniBatchSize = accumulateMiniBatch(model, miniBatchSize, miniBatch, currentMiniBatchSize, modelFeatureMapper, record);
                if (currentMiniBatchSize == miniBatchSize) {
                    // miniBatch ready, let's fit the model:
                    calibrationModel.fit(miniBatch);
                    score = calibrationModel.score();
                    currentMiniBatchSize = 0;
                    pgReadWrite.update();
                } else {
                    pgReadWrite.lightUpdate();
                }
                iter++;
            }
            pgReadWrite.stop();
            if (score < bestScore) {
                bestScore = score;
                saver.saveModel(calibrationModel, modelPrefix + "Calibrated", score);
                System.out.println("Saving best score model.. score=" + bestScore);
                lastIter = iter;

            }
            reader.close();
            reader = new RecordReader(calibrationDatafile);

        }
        reader.close();
    }

    int indices[] = {0, 0};

    /**
     * Construct a balanced mini-batch with representative of correct and incorrect predictions:
     *
     * @param model
     * @param miniBatchSize
     * @param miniBatch
     * @param currentMiniBatchSize
     * @param modelFeatureMapper
     * @param record
     * @return
     */
    private int accumulateMiniBatch(MultiLayerNetwork model, int miniBatchSize, DataSet miniBatch, int currentMiniBatchSize,
                                    FeatureMapper modelFeatureMapper, BaseInformationRecords.BaseInformation record) {
        int indexOfNewRecordInMinibatch = currentMiniBatchSize;
        ProtoPredictor predictor = new ProtoPredictor(model, modelFeatureMapper);
        boolean isTumor = record.getSamplesList().get(1).getIsTumor();
        ProtoPredictor.Prediction prediction = predictor.mutPrediction(record);
        boolean isCorrect = prediction.isCorrect(isTumor);
        // consider only errors predicting mutation (false positive predictions):
        if (!prediction.isMutated()) {

            return indexOfNewRecordInMinibatch;
        }
        int half = miniBatchSize / 2;
        if (isCorrect && numCorrectInMiniBatch(miniBatch, indexOfNewRecordInMinibatch) > half) {
            return indexOfNewRecordInMinibatch;
        }
        if (!isCorrect && numIncorrectInMiniBatch(miniBatch, indexOfNewRecordInMinibatch) > half) {
            return indexOfNewRecordInMinibatch;
        }
        // calculate the model output and use as features for the calibration model:
        FloatArrayList floats = getModelInternalActivations(model, modelFeatureMapper, record, indexOfNewRecordInMinibatch);
        indices[0] = indexOfNewRecordInMinibatch;
        for (int i = 0; i < floats.size(); i++) {
            indices[1] = i;
            miniBatch.getFeatures().putScalar(indices, floats.getFloat(i));
        }
        // run the model again to obtain the predicted label (possibility to optimize):
        // set calibration model as 1 when the prediction was correct, 1 otherwise.
        indices[1] = 0;
        miniBatch.getLabels().putScalar(indices, isCorrect ? 1 : 0);
        return ++indexOfNewRecordInMinibatch;
    }

    int[] indices2 = {0, 0};

    private int numCorrectInMiniBatch(DataSet miniBatch, int maxIndex) {
        int numCorrect = 0;
        indices2[1] = 0;
        for (int i = 0; i < maxIndex; i++) {
            indices2[0] = i;
            numCorrect += miniBatch.getLabels().getFloat(indices2) == 0 ? 1 : 0;
        }
        return numCorrect;
    }

    private int numIncorrectInMiniBatch(DataSet miniBatch, int maxIndex) {

        int numIncorrect = 0;
        indices2[1] = 0;
        for (int i = 0; i < maxIndex; i++) {
            indices2[0] = i;
            numIncorrect += miniBatch.getLabels().getFloat(indices2) == 1 ? 1 : 0;
        }
        return numIncorrect;
    }


    private FloatArrayList getModelInternalActivations(MultiLayerNetwork model, FeatureMapper modelFeatureMapper, BaseInformationRecords.BaseInformation record, int indexOfNewRecordInMinibatch) {
        INDArray inputFeatures = Nd4j.zeros(1, modelFeatureMapper.numberOfFeatures());
        modelFeatureMapper.mapFeatures(record, inputFeatures, 0);

        FloatArrayList floats = new FloatArrayList();
        model.feedForward(inputFeatures).stream().forEach(indArray -> floats.addAll(FloatArrayList.wrap(indArray.data().asFloat())));
        return floats;
    }

    private int getModelActivationNumber(MultiLayerNetwork model, FeatureMapper modelFeatureMapper) {
        int numActivations = 0;
        Layer[] layers = model.getLayers();

        INDArray inputFeatures = Nd4j.zeros(1, modelFeatureMapper.numberOfFeatures());

        int sum = model.feedForward(inputFeatures, false).stream().mapToInt(indArray ->
                indArray.data().asFloat().length).sum();
        System.out.println("Number of activations: " + sum);
        return sum;
    }

    String modelDirectory;
    String modelPrefix;
    String calibrationDatafile;

    public CalibrateModelOutput(String modelDirectory, String modelPrefix, String calibrationDatafile) {
        this.modelDirectory = modelDirectory;
        this.modelPrefix = modelPrefix;
        this.calibrationDatafile = calibrationDatafile;
    }
}
