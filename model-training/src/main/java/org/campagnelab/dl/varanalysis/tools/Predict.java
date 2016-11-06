package org.campagnelab.dl.varanalysis.tools;

import it.unimi.dsi.logging.ProgressLogger;
import org.campagnelab.dl.model.utils.mappers.FeatureMapper;
import org.campagnelab.dl.model.utils.mappers.SimpleFeatureCalculator;
import org.campagnelab.dl.model.utils.models.ModelLoader;
import org.campagnelab.dl.varanalysis.learning.iterators.BaseInformationIterator;
import org.campagnelab.dl.varanalysis.stats.AUCHelper;
import org.campagnelab.dl.varanalysis.stats.AreaUnderTheROCCurve;
import org.campagnelab.dl.varanalysis.storage.RecordReader;
import org.deeplearning4j.nn.api.Model;
import org.deeplearning4j.nn.multilayer.MultiLayerNetwork;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;

/**
 * A tool to predict with a model on a .sbi file.
 * Created by fac2003 on 10/23/16.
 */
public class Predict extends AbstractTool<PredictArguments> {

    static private Logger LOG = LoggerFactory.getLogger(Predict.class);

    public static void main(String[] args) {

        Predict predict = new Predict();
        predict.parseArguments(args, "Predict", predict.createArguments());
        predict.execute();
    }


    public PredictArguments args() {
        return (PredictArguments) arguments;
    }

    @Override
    public PredictArguments createArguments() {
        return new PredictArguments();
    }

    @Override
    public void execute() {
        try {
            File modelPath = new File(args().modelPath);
            String modelName = modelPath.getName();
            PrintWriter resultWriter = null;
            if (args().toFile) {
                String resultPath = "tests/" + modelName + "/";
                File dir = new File(resultPath);
                // attempt to create the directory here
                dir.mkdirs();

                String resultFilename = resultPath + args().modelName + "-" + args().type + ".tsv";
                System.out.println("Writing predictions to " + resultFilename);
                resultWriter = new PrintWriter(resultFilename, "UTF-8");
            } else {
                resultWriter = new PrintWriter(System.out);
            }

            printPredictions(args().modelName, args().modelPath, args().testSet, resultWriter);
        } catch (IOException e) {
            System.err.println("Error predicting in %s dataset.");
            e.printStackTrace();
            System.exit(1);
        }

    }

    private void printPredictions(String prefix, String modelPath, String evaluationDataFilename,
                                  PrintWriter resultsPath) throws IOException {


        ModelLoader modelLoader = new ModelLoader(modelPath);
        RecordReader reader = new RecordReader(evaluationDataFilename);
        FeatureMapper featureMapper = modelLoader.loadFeatureMapper(reader.getProperties());
        boolean isTrio;
        if (featureMapper.getClass().getCanonicalName().contains("Trio")) {
            //we have a trio mapper, need to output features for a third sample
            isTrio = true;
            System.out.println("setting output to trio mode");
        }
        MultiLayerNetwork model = modelLoader.loadMultiLayerNetwork(prefix);
        if (model == null) {
            System.err.println("Cannot load model with prefix: " + prefix);
            System.exit(1);
        }

        // initialize results printer
        PrintWriter results = resultsPath;
        results.append("index\ttrueLabel\tprobabilityYes\tprobabilityNo\tcorrectness").append("\n");

        ProgressLogger pgReadWrite = new ProgressLogger(LOG);
        pgReadWrite.itemsName = prefix;
        pgReadWrite.expectedUpdates = Math.min(args().scoreN, reader.getTotalRecords());
        pgReadWrite.displayFreeMemory = true;
        pgReadWrite.start();


        AreaUnderTheROCCurve aucLossCalculator = new AreaUnderTheROCCurve(args().numRecordsForAUC);
        int index = 0;
        SimpleFeatureCalculator labelMapper = new SimpleFeatureCalculator();
        BaseInformationIterator iterator = new BaseInformationIterator(evaluationDataFilename, args().miniBatchSize, featureMapper, labelMapper);
        AUCHelper helper = new AUCHelper();

        double auc =
                helper.estimate(iterator, model,
                        args().numRecordsForAUC,
                        prediction -> {
                            String correctness = (prediction.predictedLabelYes > prediction.predictedLabelNo && prediction.trueLabelYes == 1f ||
                                    prediction.predictedLabelNo > prediction.predictedLabelYes && prediction.trueLabelYes == 0f) ? "correct" : "wrong";

                            if (doOuptut(correctness, args(), Math.max(prediction.predictedLabelNo, prediction.predictedLabelYes))) {
                                results.printf("%d\t%f\t%f\t%f\t%s%n", prediction.index, prediction.trueLabelYes, prediction.predictedLabelNo,
                                        prediction.predictedLabelYes, correctness);
                            }

                            //convert true label to the convention used by auc calculator: negative true label=labelNo.
                            pgReadWrite.lightUpdate();

                        },
                /* stop if */ nProcessed -> nProcessed > args().scoreN

                );
        System.out.println("AUC on " + prefix + "=" + auc);
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