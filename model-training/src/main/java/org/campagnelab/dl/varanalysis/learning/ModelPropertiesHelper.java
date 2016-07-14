package org.campagnelab.dl.varanalysis.learning;

import org.campagnelab.dl.model.utils.mappers.FeatureMapper;
import sun.misc.ProxyGenerator;

import java.io.File;
import java.io.FileOutputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Properties;

/**
 * Model properties helper. Records properties about the model and training process in the model directory.
 * Created by fac2003 on 7/12/16.
 */
public class ModelPropertiesHelper {

    private FeatureMapper featureCalculator;
    private double learningRate;
    private int miniBatchSize;
    private int numHiddenNodes;
    private int numEpochs;
    private int seed;
    private long time;
    private int numTrainingSets;
    private double bestScore;
    private String lossFunction;

    public void writeProperties(String modelDirectory) throws IOException {
        //write properties file to model folder
        Properties modelProp = new Properties();
        modelProp.setProperty("mapper", featureCalculator.getClass().getName());
        modelProp.setProperty("learningRate", Double.toString(learningRate));
        modelProp.setProperty("time", Long.toString(time));
        modelProp.setProperty("miniBatchSize", Integer.toString(miniBatchSize));
        modelProp.setProperty("numEpochs", Integer.toString(numEpochs));
        modelProp.setProperty("numTrainingSets", Integer.toString(numTrainingSets));
        modelProp.setProperty("randSeed", Integer.toString(seed));
        modelProp.setProperty("earlyStopCondition", Integer.toString(-1));
        modelProp.setProperty("hiddenNodes", Integer.toString(numHiddenNodes));
        modelProp.setProperty("bestScore", Double.toString(bestScore));
        modelProp.setProperty("lossFunction", lossFunction);

        File file = new File(modelDirectory + "/config.properties");
        FileWriter fileWriter = new FileWriter(file);
        modelProp.store(fileWriter, "The settings used to generate this model");
    }

    public void setFeatureCalculator(FeatureMapper featureCalculator) {
        this.featureCalculator = featureCalculator;
    }

    public void setLearningRate(double learningRate) {
        this.learningRate = learningRate;
    }

    public void setNumHiddenNodes(int numHiddenNodes) {
        this.numHiddenNodes = numHiddenNodes;
    }

    public void setMiniBatchSize(int miniBatchSize) {
        this.miniBatchSize = miniBatchSize;
    }

    public void setNumEpochs(int numEpochs) {
        this.numEpochs = numEpochs;
    }

    public void setSeed(int seed) {
        this.seed = seed;
    }

    public void setTime(long time) {
        this.time = time;
    }

    public void setNumTrainingSets(int numTrainingSets) {
        this.numTrainingSets = numTrainingSets;
    }

    public void setBestScore(double bestScore) {
        this.bestScore = bestScore;
    }

    public void setLossFunction(String lossFunction) {
        this.lossFunction = lossFunction;
    }
}
