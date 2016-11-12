package org.campagnelab.dl.varanalysis.tools;

import com.google.common.collect.Iterables;
import it.unimi.dsi.logging.ProgressLogger;
import org.campagnelab.dl.model.utils.mappers.FeatureMapper;
import org.campagnelab.dl.model.utils.models.ModelLoader;
import org.campagnelab.dl.varanalysis.learning.domains.DomainDescriptor;
import org.campagnelab.dl.varanalysis.learning.PredictWithModel;
import org.campagnelab.dl.varanalysis.learning.domains.DomainDescriptorLoader;
import org.campagnelab.dl.varanalysis.learning.domains.predictions.RegressionPrediction;
import org.campagnelab.dl.varanalysis.learning.domains.predictions.SomaticFrequencyPrediction;
import org.campagnelab.dl.varanalysis.protobuf.BaseInformationRecords;
import org.campagnelab.dl.varanalysis.stats.AreaUnderTheROCCurve;
import org.campagnelab.dl.varanalysis.learning.domains.predictions.BinaryClassPrediction;
import org.campagnelab.dl.varanalysis.storage.RecordReader;
import org.deeplearning4j.nn.api.Model;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.Properties;

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
        Properties sbiProperties = reader.getProperties();
        FeatureMapper featureMapper = modelLoader.loadFeatureMapper(sbiProperties);
        boolean isTrio;
        if (featureMapper.getClass().getCanonicalName().contains("Trio")) {
            //we have a trio mapper, need to output features for a third sample
            isTrio = true;
            System.out.println("setting output to trio mode");
        }
        Model model = modelLoader.loadModel(prefix);
        if (model == null) {
            System.err.println("Cannot load model with prefix: " + prefix);
            System.exit(1);
        }

        // initialize results printer
        PrintWriter results = resultsPath;

        ProgressLogger pgReadWrite = new ProgressLogger(LOG);
        pgReadWrite.itemsName = prefix;
        pgReadWrite.expectedUpdates = Math.min(args().scoreN, reader.getTotalRecords());
        pgReadWrite.displayFreeMemory = true;
        pgReadWrite.start();
        DomainDescriptor<BaseInformationRecords.BaseInformation> domainDescriptor = DomainDescriptorLoader.load(modelPath);
        PredictWithModel<BaseInformationRecords.BaseInformation> predictor = new PredictWithModel<>(domainDescriptor);

        Iterable<BaseInformationRecords.BaseInformation> apply = domainDescriptor.getRecordIterable().apply(evaluationDataFilename);
        Iterable<BaseInformationRecords.BaseInformation> it = Iterables.limit(apply, args().numRecordsForAUC);
        AreaUnderTheROCCurve aucLossCalculator = new AreaUnderTheROCCurve(args().numRecordsForAUC);

        boolean hasSomaticFrequency = domainDescriptor.hasOutput("somaticFrequency");
        String somaticFrequencyColumns = hasSomaticFrequency ? "\ttrueSomaticFrequency\tpredictedSomaticFrequency\trmseSomaticFrequency" : "";
        results.append("index\ttrueLabel\tprobabilityYes\tprobabilityNo\tcorrectness" + somaticFrequencyColumns).append("\n");

        predictor.makePredictions(it.iterator(),
                model,
                predictionList -> {
                    // List contains at least one prediction: isSomaticMutation. It may also contain the prediction of
                    // somaticFrequency. In the second element, when the model is a computational graph with two outputs.
                    BinaryClassPrediction isSomaticMutation = (BinaryClassPrediction) predictionList.get(0);
                    String somaticFrequencyText = "";
                    if (predictionList.size() >= 2) {
                        SomaticFrequencyPrediction somaticFrequency = (SomaticFrequencyPrediction) predictionList.get(1);
                        double rmse = somaticFrequency.trueValue==null ? -1: Math.sqrt(Math.pow(somaticFrequency.trueValue - somaticFrequency.predictedValue, 2));
                        somaticFrequencyText += String.format("\t%f\t%f\t%f", somaticFrequency.trueValue, somaticFrequency.predictedValue,
                                rmse);
                    }
                    String correctness = (isSomaticMutation.predictedLabelYes > isSomaticMutation.predictedLabelNo && isSomaticMutation.trueLabelYes == 1f ||
                            isSomaticMutation.predictedLabelNo > isSomaticMutation.predictedLabelYes && isSomaticMutation.trueLabelYes == 0f) ? "correct" : "wrong";

                    if (doOuptut(correctness, args(), Math.max(isSomaticMutation.predictedLabelNo, isSomaticMutation.predictedLabelYes))) {
                        results.printf("%d\t%f\t%f\t%f\t%s%s%n", isSomaticMutation.index, isSomaticMutation.trueLabelYes, isSomaticMutation.predictedLabelYes,
                                isSomaticMutation.predictedLabelNo, correctness, somaticFrequencyText);
                    }

                    aucLossCalculator.observe(isSomaticMutation.predictedLabelYes, isSomaticMutation.trueLabelYes - 0.5);

                    //convert true label to the convention used by auc calculator: negative true label=labelNo.
                    pgReadWrite.lightUpdate();

                },
                /* stop if */ nProcessed -> nProcessed > args().scoreN
        );


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