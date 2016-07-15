package org.campagnelab.dl.varanalysis.learning;

import org.apache.commons.io.FileUtils;
import org.campagnelab.dl.model.utils.mappers.FeatureMapper;
import org.deeplearning4j.earlystopping.saver.LocalFileModelSaver;
import org.deeplearning4j.nn.conf.MultiLayerConfiguration;
import org.deeplearning4j.nn.multilayer.MultiLayerNetwork;
import org.nd4j.linalg.api.ndarray.INDArray;
import org.nd4j.linalg.factory.Nd4j;

import java.io.*;
import java.lang.reflect.Constructor;
import java.lang.reflect.InvocationTargetException;
import java.util.Properties;

/**
 * Created by fac2003 on 7/15/16.
 */
public class ModelLoader {
    String modelPath;

    public ModelLoader(String modelPath) {
        this.modelPath = modelPath;
    }

    void writeTestCount(long testRecordCount){
        try {
            FileInputStream input = new FileInputStream(modelPath + "/config.properties");
            // load a properties file
            Properties prop = new Properties();
            prop.load(input);
            input.close();
            // get the property value and print it out
            prop.setProperty("testRecordCount",Long.toString(testRecordCount));
            FileOutputStream output = new FileOutputStream(modelPath + "/config.properties");
            prop.store(output,"total testRecords added other settings, for use in statistics");
            output.close();
        } catch (FileNotFoundException e) {
            e.printStackTrace();
        } catch (IOException e) {
            e.printStackTrace();
        }


    }

    FeatureMapper loadFeatureMapper() {

        try {
            FileInputStream input = new FileInputStream(modelPath + "/config.properties");
            // load a properties file
            Properties prop = new Properties();
            prop.load(input);

            // get the property value and print it out
            String mapperName = prop.getProperty("mapper");


            ClassLoader classLoader = this.getClass().getClassLoader();

            // Load the target class using its binary name
            Class loadedMyClass = classLoader.loadClass(mapperName);

            System.out.println("Loaded class name: " + loadedMyClass.getName());

            // Create a new instance from the loaded class
            Constructor constructor = loadedMyClass.getConstructor();
            FeatureMapper featureMapper = (FeatureMapper) constructor.newInstance();
            return featureMapper;
        } catch (Exception e) {

            throw new RuntimeException("Unable to load model properties and initialize feature mapper.", e);
        }
    }

    public MultiLayerNetwork loadModel(String modelNamePrefix) throws IOException {

        MultiLayerNetwork model = null;
        if (!new File(modelPath + String.format("/%sModelParams.bin", modelNamePrefix)).isFile()) {
            model = new LocalFileModelSaver(modelPath).getBestModel();
        } else {
            //Load parameters from disk:
            INDArray newParams;
            DataInputStream dis = new DataInputStream(new FileInputStream(modelPath + String.format("/%sModelParams.bin", modelNamePrefix)));
            newParams = Nd4j.read(dis);

            //Load network configuration from disk:
            MultiLayerConfiguration confFromJson =
                    MultiLayerConfiguration.fromJson(FileUtils.readFileToString(new File(modelPath + String.format("/%sModelConf.json", modelNamePrefix))));

            //Create a MultiLayerNetwork from the saved configuration and parameters
            model = new MultiLayerNetwork(confFromJson);
            model.init();
            model.setParameters(newParams);
        }
        return model;
    }
}
