package org.campagnelab.dl.varanalysis.learning;

import it.unimi.dsi.logging.ProgressLogger;
import org.campagnelab.dl.model.utils.BayesCalibrator;
import org.campagnelab.dl.model.utils.CalcCalibrator;
import org.campagnelab.dl.model.utils.ProtoPredictor;
import org.campagnelab.dl.model.utils.mappers.FeatureMapper;
import org.campagnelab.dl.model.utils.mappers.QualityFeatures;
import org.campagnelab.dl.varanalysis.learning.calibrate.CalibratingModel;
import org.campagnelab.dl.varanalysis.protobuf.BaseInformationRecords;
import org.campagnelab.dl.varanalysis.stats.AreaUnderTheROCCurve;
import org.campagnelab.dl.varanalysis.storage.RecordReader;
import org.deeplearning4j.nn.multilayer.MultiLayerNetwork;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.io.*;
import java.util.Arrays;
import java.util.stream.IntStream;

/**
 * Use a trained neural network model to predict mutations.
 * <p>
 * Created by rct2002 on 6/8/16.
 *
 * @author Remi Torracinta
 */
public class PredictMutationsV9 extends AbstractPredictMutations {
    static private Logger LOG = LoggerFactory.getLogger(PredictMutationsV9.class);

    final static String TIME = "1469215306035";
    final boolean SKIP0COUNTS = true;
    private ModelLoader modelLoader;

    CalcCalibrator calculator;

    // String[] unmutFilenames = new String[]{"unmut_genotypes_test_proto_VN.parquet", "training_batch/genotypes_proto_" + version + "_randomized_mutated.parquet"};
    // String[] mutFilenames = new String[]{"mutated-MHFC-13-CTL_B_NK.parquet", "training_batch/genotypes_proto_" + version + "_randomized_mutated.parquet"};
    FeatureMapper featureMapper;// = new FeatureMapperV9();
    private int scoreN = Integer.MAX_VALUE;



    public PredictMutationsV9() {
    }


    public static void main(String[] args) throws Exception {
        String datasetPath;
        if (args.length == 0) {
            datasetPath = "sample_data/profiling_data/";
        } else {
            datasetPath = args[0];
        }
        double learningRate = 0.1;
        int miniBatchSize = 100;
        System.out.println("time: " + TIME);
        String attempt = "models/" + TIME;
        String modelDir = "models/" + TIME;
        PredictMutationsV9 predictor = new PredictMutationsV9();
        //   predictor.printPredictions("best");
        String type=null;
        for (String item : args) {
            if (!new File(item).exists()) {
                type=item;
                System.out.println("Will process files of type "+type);
                System.out.flush();
            } else{
                if (type==null) {
                    System.err.println("Invalid syntax, dataset type must be specified: type file+");
                    System.exit(1);
                }else{
                    datasetPath=item;
                }
                predictor.printPredictions("bestAUC", modelDir,datasetPath , "tests/" + TIME + "/", type);
            }
        }

        System.out.println(attempt);

    }


    public String featuresToString(BaseInformationRecords.BaseInformation pos) {
        //indels not handled
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

        String features = (pos.hasFrequencyOfMutation() ? pos.getFrequencyOfMutation() : "") + "\t"
                + (pos.hasMutatedBase() ? pos.getMutatedBase() : "") + "\t"
                + pos.getReferenceIndex() + "\t"
                + pos.getPosition() + "\t"
                + pos.getReferenceBase() + "\t"
                + Arrays.toString(s1Counts) + "\t"
                + Arrays.toString(s2Counts) + "\t"
                + Arrays.toString(s1Scores) + "\t"
                + Arrays.toString(s2Scores) + "\t"
                + Integer.toString(IntStream.of(s1Counts).sum() + IntStream.of(s2Counts).sum());
        return features;
    }

    public void printPredictions(String prefix, String modelPath, String evaluationDataFilename,
                                 String resultsPath, String typeOfTestFile) throws Exception {


        modelLoader = new ModelLoader(modelPath);
        featureMapper = modelLoader.loadFeatureMapper();
        calculator = new BayesCalibrator(modelPath,false);
        MultiLayerNetwork model = modelLoader.loadModel(prefix);
        MultiLayerNetwork calibratingModel = modelLoader.loadModel(prefix + "Calibrated");
        if (cmodel == null && calibratingModel != null) {
            cmodel = new CalibratingModel(model, featureMapper, calibratingModel);
        }

        File dir = new File(resultsPath);
        // attempt to create the directory here
        dir.mkdir();

        //initialize results printer

        PrintWriter results = new PrintWriter(resultsPath + prefix + "-" + typeOfTestFile, "UTF-8");
        writeHeader(results);

        //may need to adjust batch size and write outputs piecewise if test sets are very large
        //BaseInformationIterator baseIter = new BaseInformationIterator(testsetPath, Integer.MAX_VALUE, new FeatureMapperV2(), new SimpleFeatureCalculator());
        RecordReader reader = new RecordReader(evaluationDataFilename);


        //DataSet ds = baseIter.next();
//set up logger
        ProgressLogger pgReadWrite = new ProgressLogger(LOG);
        pgReadWrite.itemsName = prefix + "/" + typeOfTestFile;
        pgReadWrite.expectedUpdates = reader.getTotalRecords();
        pgReadWrite.displayFreeMemory = true;
        pgReadWrite.start();

        AreaUnderTheROCCurve aucLossCalculator = new AreaUnderTheROCCurve();
        int index = 0;

        for (BaseInformationRecords.BaseInformation record : reader) {
            //don't bother trying to make predictions when a sample has 0 counts. model outputs nan apparently.
            if (SKIP0COUNTS) {
                boolean bothHaveCount = true;
                for (BaseInformationRecords.SampleInfo sample : record.getSamplesList()){
                    boolean hasCount = false;
                    for (BaseInformationRecords.CountInfo count : sample.getCountsList()){
                        if (count.getGenotypeCountReverseStrand() > 0 || count.getGenotypeCountForwardStrand() > 0){
                            hasCount = true;
                            break;
                        }
                    }
                    if (!hasCount){
                        bothHaveCount = false;
                        break;
                    }
                }
                if (!bothHaveCount){
                    continue;
                }
            }
            writeRecordResult(model, calibratingModel, results, featureMapper, pgReadWrite, record, aucLossCalculator, calculator);
            index++;
            if (index > scoreN) break;
        }
        System.out.println("AUC on " + prefix + "=" + aucLossCalculator.evaluateStatistic());
        results.close();
        pgReadWrite.stop();

        //write sorted trees and total count to file:
        if ("test".equals(typeOfTestFile)) { //only save for test set
            //write record count to prop file
            calculator.save();
            modelLoader.writeTestCount(reader.getTotalRecords());

        }

    }

}