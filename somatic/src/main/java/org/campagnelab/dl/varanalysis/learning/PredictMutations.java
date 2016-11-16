package org.campagnelab.dl.varanalysis.learning;

import it.unimi.dsi.logging.ProgressLogger;
import org.campagnelab.dl.varanalysis.utils.BayesCalibrator;
import org.campagnelab.dl.varanalysis.utils.CalcCalibrator;
import org.campagnelab.dl.varanalysis.utils.ProtoPredictor;
import org.campagnelab.dl.varanalysis.mappers.FeatureMapper;
import org.campagnelab.dl.varanalysis.mappers.QualityFeatures;
import org.campagnelab.dl.varanalysis.models.ModelLoader;
import org.campagnelab.dl.varanalysis.learning.calibrate.CalibratingModel;
import org.campagnelab.dl.varanalysis.protobuf.BaseInformationRecords;
import org.campagnelab.dl.varanalysis.stats.AreaUnderTheROCCurve;
import org.campagnelab.dl.varanalysis.storage.RecordReader;
import org.deeplearning4j.nn.multilayer.MultiLayerNetwork;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.io.File;
import java.io.PrintWriter;
import java.util.Arrays;
import java.util.stream.IntStream;

/**
 * Use a trained neural network model to predict mutations.
 * <p>
 * Created by rct2002 on 6/8/16.
 * @deprecated  Use Predict instead.
 * @author Remi Torracinta
 */
public class PredictMutations extends AbstractPredictMutations {
    static private Logger LOG = LoggerFactory.getLogger(PredictMutations.class);


    final boolean SKIP0COUNTS = true;

    //will be adjusted if model's loaded featuremapper is for trios. don't manually change
    boolean isTrio = false;

    private ModelLoader modelLoader;

    CalcCalibrator calculator;

    // String[] unmutFilenames = new String[]{"unmut_genotypes_test_proto_VN.parquet", "training_batch/genotypes_proto_" + version + "_randomized_mutated.parquet"};
    // String[] mutFilenames = new String[]{"mutated-MHFC-13-CTL_B_NK.parquet", "training_batch/genotypes_proto_" + version + "_randomized_mutated.parquet"};
    FeatureMapper featureMapper;// = new FeatureMapperV9();


    public PredictMutations(PredictionArguments arguments) {
        super(arguments);
    }


    public static void main(String[] args) throws Exception {
        PredictionArguments arguments = parseArguments(args, "PredictMutations");
        PredictMutations predictor = new PredictMutations(arguments);
        System.out.println("modelName: " + predictor.modelName);


        predictor.printPredictions(predictor.version, arguments.modelPath, predictor.testSet, "tests/" + predictor.modelName + "/");

        System.out.println(arguments.modelPath);

    }

    //no starting tab, no ending tab
    public String featuresToString(BaseInformationRecords.BaseInformation pos, boolean longReport) {

        int[] s1Counts = new int[10];
        int[] s2Counts = new int[10];
        for (int i = 0; i < 5; i++) {
            s1Counts[i] = pos.getSamples(0).getCounts(i).getGenotypeCountForwardStrand();
            s1Counts[i + 5] = pos.getSamples(0).getCounts(i).getGenotypeCountReverseStrand();
            s2Counts[i] = pos.getSamples(1).getCounts(i).getGenotypeCountForwardStrand();
            s2Counts[i + 5] = pos.getSamples(1).getCounts(i).getGenotypeCountReverseStrand();
        }
        float[] s1Scores = new float[10];
        float[] s2Scores = new float[10];
        for (int i = 0; i < 5; i++) {
            s1Scores[i] = QualityFeatures.avgQuality(ProtoPredictor.expandFreq(pos.getSamples(0).getCounts(i).getQualityScoresForwardStrandList()));
            s1Scores[i + 5] = QualityFeatures.avgQuality(ProtoPredictor.expandFreq(pos.getSamples(0).getCounts(i).getQualityScoresReverseStrandList()));
            s2Scores[i] = QualityFeatures.avgQuality(ProtoPredictor.expandFreq(pos.getSamples(1).getCounts(i).getQualityScoresForwardStrandList()));
            s2Scores[i + 5] = QualityFeatures.avgQuality(ProtoPredictor.expandFreq(pos.getSamples(1).getCounts(i).getQualityScoresReverseStrandList()));
        }
        String s3CountsString = "";
        String s3ScoresString = "";
        int s3CountsSum = 0;
        if (isTrio) {
            int[] s3Counts = new int[10];
            for (int i = 0; i < 5; i++) {
                s3Counts[i] = pos.getSamples(2).getCounts(i).getGenotypeCountForwardStrand();
                s3Counts[i + 5] = pos.getSamples(2).getCounts(i).getGenotypeCountReverseStrand();
            }
            float[] s3Scores = new float[10];
            for (int i = 0; i < 5; i++) {
                s3Scores[i] = QualityFeatures.avgQuality(ProtoPredictor.expandFreq(pos.getSamples(2).getCounts(i).getQualityScoresForwardStrandList()));
                s3Scores[i + 5] = QualityFeatures.avgQuality(ProtoPredictor.expandFreq(pos.getSamples(2).getCounts(i).getQualityScoresReverseStrandList()));
            }
            s3CountsString = "\t" + Arrays.toString(s3Counts).replaceAll(" ", "");
            s3ScoresString = "\t" + Arrays.toString(s3Scores).replaceAll(" ", "");
            s3CountsSum = IntStream.of(s3Counts).sum();
        }
        String features;
        String refId=(pos.hasReferenceId()?pos.getReferenceId():Integer.toString(pos.getReferenceIndex()));
        if (!longReport) {
            features = (pos.hasFrequencyOfMutation() ? pos.getFrequencyOfMutation() : "") + "\t"
                    + refId + "\t"
                    + pos.getPosition() + "\t"
                    + pos.getReferenceBase() + "\t"
                    + Arrays.toString(s1Scores).replaceAll(" ", "") + "\t"
                    + Arrays.toString(s2Scores).replaceAll(" ", "")
                    + s3ScoresString + "\t"
                    + Integer.toString(IntStream.of(s1Counts).sum() + IntStream.of(s2Counts).sum() + s3CountsSum);

        } else {
            features = (pos.hasFrequencyOfMutation() ? pos.getFrequencyOfMutation() : "") + "\t"
                    + refId + "\t"
                    + pos.getPosition() + "\t"
                    + pos.getReferenceBase() + "\t"
                    + Arrays.toString(s1Scores).replaceAll(" ", "") + "\t"
                    + Arrays.toString(s2Scores).replaceAll(" ", "")
                    + s3ScoresString + "\t"
                    + Integer.toString(IntStream.of(s1Counts).sum() + IntStream.of(s2Counts).sum() + s3CountsSum)
                    + (pos.hasMutatedBase() ? pos.getMutatedBase() : "") + "\t"
                    + Arrays.toString(s1Counts).replaceAll(" ", "") + "\t"
                    + Arrays.toString(s2Counts).replaceAll(" ", "")
                    + s3CountsString;

        }
        return features;
    }

    public void printPredictions(String prefix, String modelPath, String evaluationDataFilename,
                                 String resultsPath) throws Exception {


        modelLoader = new ModelLoader(modelPath);
        RecordReader reader = new RecordReader(evaluationDataFilename);
        featureMapper = modelLoader.loadFeatureMapper(reader.getProperties());
        if (featureMapper.getClass().getCanonicalName().contains("Trio")) {
            //we have a trio mapper, need to output features for a third sample
            isTrio = true;
            System.out.println("setting output to trio mode");
        }
        calculator = new BayesCalibrator(modelPath, prefix, false);
        MultiLayerNetwork model = modelLoader.loadMultiLayerNetwork(prefix);
        if (model == null) {
            System.err.println("Cannot load model with prefix: " + prefix);
            System.exit(1);
        }

        MultiLayerNetwork calibratingModel = modelLoader.loadMultiLayerNetwork(prefix + "Calibrated");
        if (cmodel == null && calibratingModel != null) {
            cmodel = new CalibratingModel(model, featureMapper, calibratingModel);
        }

        File dir = new File(resultsPath);
        // attempt to create the directory here
        dir.mkdir();

        //initialize results printer

        String resultFilename = resultsPath + prefix + "-" + type + ".tsv";
        PrintWriter results = new PrintWriter(resultFilename, "UTF-8");
        writeHeader(results, isTrio);

        System.out.println("Writing predictions to " + resultFilename);
        //may need to adjust batch size and write outputs piecewise if test sets are very large
        //BaseInformationIterator baseIter = new BaseInformationIterator(testsetPath, Integer.MAX_VALUE, new FeatureMapperV2(), new SimpleFeatureCalculator());


        //DataSet ds = baseIter.next();
//set up logger
        ProgressLogger pgReadWrite = new ProgressLogger(LOG);
        pgReadWrite.itemsName = prefix;
        pgReadWrite.expectedUpdates = Math.min(arguments.scoreN,reader.getTotalRecords());
        pgReadWrite.displayFreeMemory = true;
        pgReadWrite.start();

        AreaUnderTheROCCurve aucLossCalculator = new AreaUnderTheROCCurve(50000);
        int index = 0;

        for (BaseInformationRecords.BaseInformation record : reader) {
            //don't bother trying to make predictions when a sample has 0 counts. model outputs nan apparently.
            if (SKIP0COUNTS) {
                boolean bothHaveCount = true;
                for (BaseInformationRecords.SampleInfo sample : record.getSamplesList()) {
                    boolean hasCount = false;
                    for (BaseInformationRecords.CountInfo count : sample.getCountsList()) {
                        if (count.getGenotypeCountReverseStrand() > 0 || count.getGenotypeCountForwardStrand() > 0) {
                            hasCount = true;
                            break;
                        }
                    }
                    if (!hasCount) {
                        bothHaveCount = false;
                        break;
                    }
                }
                if (!bothHaveCount) {
                    continue;
                }
            }
            writeRecordResult(model, calibratingModel, results, featureMapper, pgReadWrite, record, aucLossCalculator, calculator, isTrio);
            index++;
            if (index > scoreN) break;
        }
        System.out.println("AUC on " + prefix + "=" + aucLossCalculator.evaluateStatistic());
        results.close();
        pgReadWrite.stop();


        calculator.save();
        modelLoader.writeTestCount(reader.getTotalRecords());


    }

}