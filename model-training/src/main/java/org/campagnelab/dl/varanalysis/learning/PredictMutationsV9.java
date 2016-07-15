package org.campagnelab.dl.varanalysis.learning;

import it.unimi.dsi.logging.ProgressLogger;
import org.apache.commons.io.FileUtils;
import org.campagnelab.dl.model.utils.ProtoPredictor;
import org.campagnelab.dl.model.utils.mappers.FeatureMapper;
import org.campagnelab.dl.model.utils.mappers.QualityFeatures;
import org.campagnelab.dl.varanalysis.protobuf.BaseInformationRecords;
import org.campagnelab.dl.varanalysis.stats.AreaUnderTheROCCurve;
import org.campagnelab.dl.varanalysis.storage.RecordReader;
import org.deeplearning4j.earlystopping.saver.LocalFileModelSaver;
import org.deeplearning4j.nn.conf.MultiLayerConfiguration;
import org.deeplearning4j.nn.multilayer.MultiLayerNetwork;
import org.nd4j.linalg.api.ndarray.INDArray;
import org.nd4j.linalg.factory.Nd4j;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.io.*;
import java.lang.reflect.Constructor;
import java.util.Arrays;
import java.util.Properties;
import java.util.stream.IntStream;

/**
 * Train a neural network to predict mutations.
 * <p>
 * Created by rct2002 on 6/8/16.
 *
 * @author Remi Torracinta
 */
public class PredictMutationsV9 extends AbstractPredictMutations {
    static private Logger LOG = LoggerFactory.getLogger(PredictMutationsV9.class);


    final static String TIME = "1468533491636";
    private final ModelLoader modelLoader;
    String modelPath;
    String dataDirPath;
    String resultsPath;
    String version = "VN";
    String[] dataFilenames;
    boolean TEST_UNMUT = false;
    String[] unmutFilenames = new String[]{"unmut_genotypes_test_proto_VN.parquet", "training_batch/genotypes_proto_" + version + "_randomized_mutated.parquet"};
    String[] mutFilenames = new String[]{"mutated-MHFC-13-CTL_B_NK.parquet", "training_batch/genotypes_proto_" + version + "_randomized_mutated.parquet"};
    String[] resultsFileNames = new String[]{"test", "training"};
    FeatureMapper featureMapper;// = new FeatureMapperV9();
    private int scoreN = Integer.MAX_VALUE;


    public PredictMutationsV9(String modelPath, String dataDirPath, String resultsPath) {
        this.modelPath = modelPath;
        this.dataDirPath = dataDirPath;
        this.resultsPath = resultsPath;
        this.dataFilenames = TEST_UNMUT ? unmutFilenames : mutFilenames;


        modelLoader = new ModelLoader(modelPath);
        featureMapper = modelLoader.loadFeatureMapper();
    }


    public static void main(String[] args) throws IOException {
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
        PredictMutationsV9 predictor = new PredictMutationsV9(modelDir, datasetPath, "tests/" + TIME + "/");
        predictor.printPredictions("best");
        predictor.printPredictions("latest");
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

    public void printPredictions(String prefix) {
        try {
            MultiLayerNetwork model=modelLoader.loadModel(prefix);
            File dir = new File(resultsPath);
            // attempt to create the directory here
            dir.mkdir();


            for (int i = 0; i < 2; i++) {
                //initialize results printer
                String typeOfTestFile = resultsFileNames[i];
                PrintWriter results = new PrintWriter(resultsPath + prefix + "-" + typeOfTestFile, "UTF-8");
                writeHeader(results);

                //may need to adjust batch size and write outputs piecewise if test sets are very large
                //BaseInformationIterator baseIter = new BaseInformationIterator(testsetPath, Integer.MAX_VALUE, new FeatureMapperV2(), new SimpleFeatureCalculator());
                RecordReader reader = new RecordReader(dataDirPath + dataFilenames[i]);
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
                    writeRecordResult(model, results, featureMapper, pgReadWrite, record, aucLossCalculator);
                    index++;
                    if (index > scoreN) break;
                }
                System.out.println("AUC on " + prefix + "=" + aucLossCalculator.evaluateStatistic());
                results.close();
                pgReadWrite.stop();
            }


        } catch (Exception e) {
            throw new RuntimeException((e));
        }


    }

}