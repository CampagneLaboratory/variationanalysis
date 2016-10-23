package org.campagnelab.dl.varanalysis.tools;

import it.unimi.dsi.logging.ProgressLogger;
import org.campagnelab.dl.model.utils.mappers.FeatureMapper;
import org.campagnelab.dl.model.utils.mappers.SimpleFeatureCalculator;
import org.campagnelab.dl.model.utils.models.ModelLoader;
import org.campagnelab.dl.varanalysis.learning.iterators.BaseInformationIterator;
import org.campagnelab.dl.varanalysis.stats.AreaUnderTheROCCurve;
import org.campagnelab.dl.varanalysis.storage.RecordReader;
import org.deeplearning4j.nn.multilayer.MultiLayerNetwork;
import org.nd4j.linalg.api.ndarray.INDArray;
import org.nd4j.linalg.dataset.DataSet;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;

/**
 * A tool to predict with a model on a .sbi file.
 * Created by fac2003 on 10/23/16.
 */
public class Predict extends AbstractTool {

    static private Logger LOG = LoggerFactory.getLogger(Predict.class);

    public static void main(String[] args) {
        PredictArguments arguments = new PredictArguments();
        Predict predict = new Predict();
        predict.parseArguments(args, "Predict", arguments);
        predict.execute();
    }

    public PredictArguments args() {
        return (PredictArguments) arguments;
    }

    @Override
    public void execute() {
        try {
            File modelPath = new File(args().modelPath);
            String modelName = modelPath.getName();

            printPredictions(args().modelName, args().modelPath, args().testSet, "tests/" + modelName + "/");
        } catch (IOException e) {
            System.err.println("Error predicting in %s dataset.");
            e.printStackTrace();
            System.exit(1);
        }

    }

    private void printPredictions(String prefix, String modelPath, String evaluationDataFilename,
                                  String resultsPath) throws IOException {


        ModelLoader modelLoader = new ModelLoader(modelPath);
        RecordReader reader = new RecordReader(evaluationDataFilename);
        FeatureMapper featureMapper = modelLoader.loadFeatureMapper(reader.getProperties());
        boolean isTrio;
        if (featureMapper.getClass().getCanonicalName().contains("Trio")) {
            //we have a trio mapper, need to output features for a third sample
            isTrio = true;
            System.out.println("setting output to trio mode");
        }
        MultiLayerNetwork model = modelLoader.loadModel(prefix);
        if (model == null) {
            System.err.println("Cannot load model with prefix: " + prefix);
            System.exit(1);
        }
        File dir = new File(resultsPath);

        // attempt to create the directory here
        dir.mkdirs();

        //initialize results printer

        String resultFilename = resultsPath + prefix + "-" + args().type + ".tsv";
        PrintWriter results = new PrintWriter(resultFilename, "UTF-8");
        results.append("index\ttrueLabel\tprobabilityYes\tprobabilityNo\tcorrectness").append("\n");


        System.out.println("Writing predictions to " + resultFilename);

        ProgressLogger pgReadWrite = new ProgressLogger(LOG);
        pgReadWrite.itemsName = prefix;
        pgReadWrite.expectedUpdates = Math.min(args().scoreN, reader.getTotalRecords());
        pgReadWrite.displayFreeMemory = true;
        pgReadWrite.start();

        AreaUnderTheROCCurve aucLossCalculator = new AreaUnderTheROCCurve(args().numRecordsForAUC);
        int index = 0;
        SimpleFeatureCalculator labelMapper = new SimpleFeatureCalculator();
        BaseInformationIterator iterator = new BaseInformationIterator(evaluationDataFilename, args().miniBatchSize, featureMapper, labelMapper);
        int nProcessed = 0;
        while (iterator.hasNext()) {
            DataSet next = iterator.next();
            INDArray outputs = model.output(next.getFeatures());
            for (int predictionIndex = 0; predictionIndex < next.numExamples(); predictionIndex++) {
                INDArray trueLabels = next.getLabels();

                double trueLabelYes = trueLabels.getDouble(predictionIndex, 1);
                double predictedLabelNo = outputs.getDouble(predictionIndex, 0);
                double predictedLabelYes = outputs.getDouble(predictionIndex, 1);
                String correctness = (predictedLabelYes > predictedLabelNo && trueLabelYes == 1f ||
                        predictedLabelNo > predictedLabelYes && trueLabelYes == 0) ? "correct" : "wrong";
                if (doOuptut(correctness, args(), Math.max(predictedLabelNo, predictedLabelYes))) {
                    results.printf("%d\t%f\t%f\t%f\t%s%n", index++, trueLabelYes, predictedLabelNo,
                            predictedLabelYes, correctness);
                }
                //convert true label to the convention used by auc calculator: negative true label=labelNo.
                aucLossCalculator.observe(predictedLabelYes, trueLabelYes-0.5);
                //     writeRecordResultFast(model, results, featureMapper, pgReadWrite, record, aucLossCalculator, calculator, isTrio);
                pgReadWrite.lightUpdate();
            }
            nProcessed += next.numExamples();
            if (nProcessed > args().scoreN) {
                break;
            }

        }

        System.out.println("AUC on " + prefix + "=" + aucLossCalculator.evaluateStatistic());
        results.close();
        pgReadWrite.stop();

        modelLoader.writeTestCount(reader.getTotalRecords());


    }

    /**
     * Apply filters and decide if a prediction should be written to the output.
     *
     * @param correctness
     * @param args
     * @param pMax
     * @return
     */
    private boolean doOuptut(String correctness, PredictArguments args, double pMax) {
        if (args.correctnessFilter != null) {
            if (!correctness.equals(args.correctnessFilter)) {
                return false;
            }
        }
        if (pMax < args().pFilterMinimum || pMax > args().pFilterMaximum) {
            return false;
        }
        return true;
    }


}