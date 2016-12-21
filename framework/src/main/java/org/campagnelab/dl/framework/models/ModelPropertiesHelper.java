package org.campagnelab.dl.framework.models;

import org.campagnelab.dl.framework.gpu.ParameterPrecision;
import org.campagnelab.dl.framework.mappers.FeatureMapper;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Properties;

/**
 * Model properties helper. Records properties about the model and training process in the model directory.
 * Created by fac2003 on 7/12/16.
 */
public class ModelPropertiesHelper {
    Properties modelProp = new Properties();

    public ModelPropertiesHelper() {

    }

    private int numHiddenNodes;
    private String lossFunction;

    public void setPrecision(ParameterPrecision p) {
        modelProp.setProperty("precision", p.name());
    }

    public void writeProperties(String modelDirectory) throws IOException {
        //write properties file to model folder
        modelProp.setProperty("hiddenNodes", Integer.toString(numHiddenNodes));
        if (lossFunction != null) {
            modelProp.setProperty("lossFunction", lossFunction);
        }
        File file = new File(modelDirectory + "/config.properties");
        FileWriter fileWriter = new FileWriter(file);
        modelProp.store(fileWriter, "The settings used to generate this model");
    }

    public void setFeatureCalculator(FeatureMapper featureCalculator) {
        modelProp.setProperty("mapper", featureCalculator.getClass().getName());
    }

    public void setLearningRate(double learningRate) {
        modelProp.setProperty("learningRate", Double.toString(learningRate));
    }

    public void setDropoutRate(Double dropoutRate) {
        if (dropoutRate != null) {
            modelProp.setProperty("dropoutRate", Double.toString(dropoutRate));
        }
    }

    public void setRegularization(Double regularization) {
      if (regularization!=null) {
          modelProp.setProperty("regularizationRate", Double.toString(regularization));
      }
    }

    public void setNumHiddenNodes(int numHiddenNodes) {
        this.numHiddenNodes = numHiddenNodes;
    }

    public void setMiniBatchSize(int miniBatchSize) {
        modelProp.setProperty("miniBatchSize", Integer.toString(miniBatchSize));

    }

    public void setNumEpochs(int numEpochs) {
        modelProp.setProperty("numEpochs", Integer.toString(numEpochs));
    }

    public void setSeed(long seed) {
        modelProp.setProperty("randomSeed", Long.toString(seed));
    }

    public void setTime(long time) {
        modelProp.setProperty("time", Long.toString(time));
    }

    public void setNumTrainingSets(int numTrainingSets) {
        modelProp.setProperty("numTrainingSets", Integer.toString(numTrainingSets));
    }

    public void setEarlyStopCriterion(int value) {
        modelProp.setProperty("earlyStopCondition", Integer.toString(value));
    }

    public void setLossFunction(String lossFunction) {
        this.lossFunction = lossFunction;
    }

    public void put(String key, int value) {
        modelProp.setProperty(key, Integer.toString(value));
    }

    public void put(String key, String value) {
        modelProp.setProperty(key, value);
    }

    public void put(String key, boolean value) {
        modelProp.setProperty(key, Boolean.toString(value));
    }


    public void addProperties(Properties readerProperties) {
        this.modelProp.putAll(readerProperties);
    }

    public Properties getProperties() {
        return modelProp;
    }
}
