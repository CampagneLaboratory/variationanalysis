package org.campagnelab.dl.varanalysis.learning;

import com.google.protobuf.TextFormat;
import org.apache.commons.io.FileUtils;
import org.campagnelab.dl.varanalysis.learning.iterators.BaseInformationIterator;
import org.campagnelab.dl.varanalysis.learning.mappers.FeatureMapperV2;
import org.campagnelab.dl.varanalysis.learning.iterators.SimpleFeatureCalculator;
import org.campagnelab.dl.varanalysis.protobuf.BaseInformationRecords;
import org.campagnelab.dl.varanalysis.storage.RecordReader;
import org.deeplearning4j.nn.conf.MultiLayerConfiguration;
import org.deeplearning4j.nn.multilayer.MultiLayerNetwork;
import org.nd4j.linalg.api.ndarray.INDArray;
import org.nd4j.linalg.dataset.DataSet;
import org.nd4j.linalg.factory.Nd4j;

import java.util.*;

import java.io.*;

/**
 * Train a neural network to predict mutations.
 * <p>
 * Created by rct2002 on 5/23/16.
 *
 * @author Remi Torracinta
 */
public class PredictMutations {
    String modelPath;
    String dataDirPath;
    String resultsPath;
    String[] dataFilenames = new String[]{"genotypes_proto.parquet","genotypes_proto_test.parquet"};
    String[] resultsFileNames = new String[]{"training","test"};
    final String header = "mutatedLabel\tProbability\tcorrectness\tfrequency\tmutatedBase\trefIdx\tposition\treferenceBase\tsample1Counts\tsample2Counts\tsample1Scores\tsample2Scores\n";




    public PredictMutations(String modelPath, String dataDirPath, String resultsPath){
        this.modelPath = modelPath;
        this.dataDirPath = dataDirPath;
        this.resultsPath = resultsPath;

    }


    public static void main(String[] args) throws IOException {
        double learningRate = 0.05;
        int miniBatchSize = 100;
        String time = "1464722138875";
        String attempt = "batch=" + miniBatchSize + "learningRate=" + learningRate + "-time:" + time;

        PredictMutations predictor = new PredictMutations(attempt, "sample_data/protobuf/", "tests/" + time + "/");
        predictor.PrintPredictions();
    }


    String featuresToString(BaseInformationRecords.BaseInformation pos){
        //indels not handled
        int[] s1Counts = new int[10];
        int[] s2Counts = new int[10];
        for (int i = 0; i < 5; i++){
            s1Counts[i] = pos.getSamples(0).getCounts(i).getGenotypeCountForwardStrand();
            s1Counts[i+5] = pos.getSamples(0).getCounts(i).getGenotypeCountReverseStrand();
            s2Counts[i] = pos.getSamples(1).getCounts(i).getGenotypeCountForwardStrand();
            s2Counts[i+5] = pos.getSamples(1).getCounts(i).getGenotypeCountReverseStrand();
        }
        float[] s1Scores = new float[10];
        float[] s2Scores = new float[10];
        for (int i = 0; i < 5; i++){
            s1Scores[i] = pos.getSamples(0).getCounts(i).getQualityScoreForwardStrand();
            s1Scores[i+5] = pos.getSamples(0).getCounts(i).getQualityScoreReverseStrand();
            s2Scores[i] = pos.getSamples(1).getCounts(i).getQualityScoreForwardStrand();
            s2Scores[i+5] = pos.getSamples(1).getCounts(i).getQualityScoreReverseStrand();
        }

        String features = (pos.hasFrequencyOfMutation()?pos.getFrequencyOfMutation():"") + "\t"
                        + (pos.hasMutatedBase()?pos.getMutatedBase():"") + "\t"
                        + pos.getReferenceIndex() + "\t"
                        + pos.getPosition() + "\t"
                        + pos.getReferenceBase() + "\t"
                        + Arrays.toString(s1Counts) + "\t"
                        + Arrays.toString(s2Counts) + "\t"
                        + Arrays.toString(s1Scores) + "\t"
                        + Arrays.toString(s2Scores);
        return features;
    }

    public void PrintPredictions(){
        try {
            //Load parameters from disk:
            INDArray newParams;
            DataInputStream dis = new DataInputStream(new FileInputStream(modelPath + "/bestModelParams.bin"));
            newParams = Nd4j.read(dis);

            //Load network configuration from disk:
            MultiLayerConfiguration confFromJson =
                    MultiLayerConfiguration.fromJson(FileUtils.readFileToString(new File(modelPath + "/bestModelConf.json")));

            //Create a MultiLayerNetwork from the saved configuration and parameters
            MultiLayerNetwork model = new MultiLayerNetwork(confFromJson);
            model.init();
            model.setParameters(newParams);

            File dir = new File(resultsPath);
            // attempt to create the directory here
            dir.mkdir();


            for (int i = 0; i < 2; i++){
                //initialize results printer
                PrintWriter results = new PrintWriter(resultsPath+resultsFileNames[i], "UTF-8");
                results.append(header);

                //may need to adjust batch size and write outputs piecewise if test sets are very large
                //BaseInformationIterator baseIter = new BaseInformationIterator(testsetPath, Integer.MAX_VALUE, new FeatureMapperV2(), new SimpleFeatureCalculator());
                FeatureMapperV2 featureMapper = new FeatureMapperV2();
                RecordReader reader = new RecordReader(dataDirPath+dataFilenames[i]);
                //DataSet ds = baseIter.next();


                for (BaseInformationRecords.BaseInformation record : reader){
                    INDArray testFeatures = Nd4j.zeros(1, featureMapper.numberOfFeatures());
                    featureMapper.mapFeatures(record,testFeatures,0);
                    INDArray testPredicted = model.output(testFeatures,false);
                    String features = featuresToString(record);
                    //boolean
                    boolean mutated = record.getMutated();
                    float[] probabilities = testPredicted.getRow(0).data().asFloat();
                    boolean prediction = probabilities[0] > 0.5;
                    String correctness = (prediction == mutated) ? "right" : "wrong";
                    results.append((mutated?"1":"0") + "\t" + Float.toString(probabilities[0]) + "\t" + correctness + "\t" + features + "\n");
                }
                results.close();
            }


        } catch (Exception e) {
            throw new RuntimeException((e));
        }




    }

}