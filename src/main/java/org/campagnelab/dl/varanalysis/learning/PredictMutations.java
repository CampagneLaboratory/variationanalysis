package org.campagnelab.dl.varanalysis.learning;

import it.unimi.dsi.logging.ProgressLogger;
import org.apache.commons.io.FileUtils;
import org.campagnelab.dl.varanalysis.format.PosRecord;
import org.campagnelab.dl.varanalysis.format.SampleRecord;
import org.campagnelab.dl.varanalysis.learning.iterators.BaseInformationIterator;
import org.campagnelab.dl.varanalysis.learning.iterators.SimpleFeatureCalculator;
import org.campagnelab.dl.varanalysis.storage.AvroVariationParquetReader;
import org.canova.api.records.reader.RecordReader;
import org.canova.api.records.reader.impl.CSVRecordReader;
import org.canova.api.split.FileSplit;
import org.deeplearning4j.datasets.canova.RecordReaderDataSetIterator;
import org.deeplearning4j.datasets.iterator.DataSetIterator;
import org.deeplearning4j.earlystopping.saver.LocalFileModelSaver;
import org.deeplearning4j.nn.api.Layer;
import org.deeplearning4j.nn.api.OptimizationAlgorithm;
import org.deeplearning4j.nn.conf.MultiLayerConfiguration;
import org.deeplearning4j.nn.conf.NeuralNetConfiguration;
import org.deeplearning4j.nn.conf.Updater;
import org.deeplearning4j.nn.conf.layers.DenseLayer;
import org.deeplearning4j.nn.conf.layers.OutputLayer;
import org.deeplearning4j.nn.multilayer.MultiLayerNetwork;
import org.deeplearning4j.nn.weights.WeightInit;
import org.deeplearning4j.optimize.api.IterationListener;
import org.deeplearning4j.optimize.listeners.ScoreIterationListener;
import org.deeplearning4j.ui.weights.HistogramIterationListener;
import org.nd4j.linalg.api.ndarray.INDArray;
import org.nd4j.linalg.dataset.DataSet;
import org.nd4j.linalg.factory.Nd4j;
import org.nd4j.linalg.lossfunctions.LossFunctions;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.util.*;

import java.io.*;
import java.nio.charset.Charset;
import java.nio.file.Files;
import java.nio.file.Paths;

/**
 * Train a neural network to predict mutations.
 * <p>
 * Created by rct2002 on 5/23/16.
 *
 * @author Remi Torracinta
 */
public class PredictMutations {
    String modelPath;
    String testsetPath;
    String resultsPath;
    public PredictMutations(String modelPath, String testsetPath, String resultsPath){
        this.modelPath = modelPath;
        this.testsetPath = testsetPath;
        this.resultsPath = resultsPath;

    }


    public static void main(String[] args) throws IOException {
        double learningRate = 0.001;
        int miniBatchSize = 1000;
        String attempt = "batch=" + miniBatchSize + "learningRate=" + learningRate;

        PredictMutations predictor = new PredictMutations(attempt, "sample_data/genotypes_testset_randomized.parquet", "tests/results");
        predictor.PrintPredictions();
    }


    String featuresToString(PosRecord pos){
        StringBuilder sb = new StringBuilder();

        List<SampleRecord> srecs = pos.getSamples();
        for (SampleRecord srec : srecs) {
            sb.append(srec.getCounts().toString());
        }

        return sb.toString();
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



            //initialize results printer
            PrintWriter results = new PrintWriter(resultsPath, "UTF-8");

            //may need to adjust batch size and write outputs piecewise if test sets are very large
            BaseInformationIterator baseIter = new BaseInformationIterator(testsetPath, Integer.MAX_VALUE, new SimpleFeatureCalculator());
            AvroVariationParquetReader avroReader = new AvroVariationParquetReader(testsetPath);
            DataSet ds = baseIter.next();
            INDArray testFeatures = ds.getFeatures();
            INDArray testPredicted = model.output(testFeatures);
            INDArray testActual = ds.getLabels();

            //iterate over features + labels and print to tab delimited file
            for (int i = 0; i < testPredicted.rows(); i++){
                PosRecord pos = avroReader.read();
                String features = featuresToString(pos);
                //boolean
                boolean mutated = pos.getMutated();
                float[] probabilities = testPredicted.getRow(i).data().asFloat();
                //float[] labels = testActual.getRow(i).data().asFloat();
                boolean prediction = probabilities[0] > 0.5;
                String correctness = (prediction == mutated) ? "right" : "wrong";
                results.append((mutated?"1":"0") + "\t" + Float.toString(probabilities[0]) + "\t" + correctness + "\t" + features + "\n");
            }
            results.close();
        } catch (Exception e) {
            throw new RuntimeException((e));
        }




    }

}
