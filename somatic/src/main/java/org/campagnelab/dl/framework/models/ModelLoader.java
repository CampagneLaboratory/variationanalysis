package org.campagnelab.dl.framework.models;

import org.apache.commons.io.FileUtils;
import org.campagnelab.dl.framework.mappers.ConfigurableFeatureMapper;
import org.campagnelab.dl.framework.mappers.FeatureMapper;
import org.deeplearning4j.nn.api.Model;
import org.deeplearning4j.nn.conf.MultiLayerConfiguration;
import org.deeplearning4j.nn.multilayer.MultiLayerNetwork;
import org.deeplearning4j.util.ModelSerializer;
import org.nd4j.linalg.api.buffer.DataBuffer;
import org.nd4j.linalg.api.buffer.util.DataTypeUtil;
import org.nd4j.linalg.api.ndarray.INDArray;
import org.nd4j.linalg.factory.Nd4j;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.io.*;
import java.lang.reflect.Constructor;
import java.nio.charset.Charset;
import java.util.Properties;

/**
 * Created by fac2003 on 7/15/16.
 */
public class ModelLoader {
    private final Properties modelProperties;
    String modelPath;
    static private Logger LOG = LoggerFactory.getLogger(ModelLoader.class);

    public ModelLoader(String modelPath) {
        this.modelPath = modelPath;

        try {
            modelProperties = loadModelProperties();
        } catch (IOException e) {
            throw new RuntimeException("Unable to load model properties at path " + modelPath);
        }
    }

    /**
     * Write the number of records used in test set to model config.properties.
     *
     * @param testRecordCount
     */
    public void writeTestCount(long testRecordCount) {
        try {
            FileInputStream input = new FileInputStream(modelPath + "/config.properties");
            // load a properties file
            Properties prop = new Properties();
            prop.load(input);
            input.close();
            // get the property value and print it out
            prop.setProperty("testRecordCount", Long.toString(testRecordCount));
            FileOutputStream output = new FileOutputStream(modelPath + "/config.properties");
            prop.store(output, "total testRecords added other settings, for use in statistics");
            output.close();
        } catch (FileNotFoundException e) {
            e.printStackTrace();
        } catch (IOException e) {
            e.printStackTrace();
        }


    }

    public FeatureMapper loadFeatureMapper(Properties sbiProperties) {

        try {

            // get the property value and print it out
            String mapperName = modelProperties.getProperty("mapper");
            ClassLoader classLoader = this.getClass().getClassLoader();

            // Load the target class using its binary name
            Class loadedMyClass = classLoader.loadClass(mapperName);

            System.out.println("Loaded class name: " + loadedMyClass.getName());

            // Create a new instance from the loaded class
            Constructor constructor = loadedMyClass.getConstructor();
            FeatureMapper featureMapper = (FeatureMapper) constructor.newInstance();
            if (featureMapper instanceof ConfigurableFeatureMapper) {
                ConfigurableFeatureMapper cfm = (ConfigurableFeatureMapper) featureMapper;
                cfm.configure(sbiProperties);
            }
            return featureMapper;
        } catch (Exception e) {

            throw new RuntimeException("Unable to load model properties and initialize feature mapper.", e);
        }
    }

    private Properties loadModelProperties() throws IOException {
        FileInputStream input = new FileInputStream(modelPath + "/config.properties");
        // load a properties file
        Properties prop = new Properties();
        prop.load(input);
        if (prop.getProperty("precision") != null && prop.getProperty("precision").equals("FP16")) {
            LOG.info("Model uses FP16 precision. Activating support.");
            DataTypeUtil.setDTypeForContext(DataBuffer.Type.HALF);
        }
        return prop;
    }

    public MultiLayerNetwork loadMultiLayerNetwork(String modelNamePrefix) throws IOException {
        Model m = loadModel(modelNamePrefix);
        if (m instanceof MultiLayerNetwork) {
            return (MultiLayerNetwork) m;
        }
        return null;
    }

    /**
     * Strip the model path from the suffix to path/modelLabel
     *
     * @param fullPath
     * @return model prefix/label or null if the full path did not point to a model file.
     */
    public static String getModelLabel(String fullPath) {
        String[] modelPathSplit = fullPath.split("/");
        String modelFileName = modelPathSplit[modelPathSplit.length - 1];
        String suffixes[] = {"-ComputationGraph.bin", "Model.bin"};
        if (modelFileName.endsWith(".bin")) {
            for (String suffix : suffixes) {
                if (modelFileName.endsWith(suffix)) {
                    return modelFileName.substring(0, modelFileName.length() - suffix.length());
                }
            }
        }
        return null;
    }

    public static String getModelPath(String fullPath) {
        if (fullPath.endsWith(".bin")) {
            return new File(fullPath).getParent();
        } else {
            return fullPath;
        }
    }

    public Model loadModel(String modelNamePrefix) throws IOException {

        Model model = null;

        String pathname = getPath(modelNamePrefix, "/%sModel.bin");
        if (new File(pathname).exists()) {
            model = ModelSerializer.restoreMultiLayerNetwork(pathname);
            return model;
        } else {
            pathname = getPath(modelNamePrefix, "/%s-ComputationGraph.bin");
            if (new File(pathname).exists()) {
                model = ModelSerializer.restoreComputationGraph(pathname);
                return model;
            }
        }

        MultiLayerNetwork net = null;

        if (!(new File(pathname).exists() || new File(getPath(modelNamePrefix, "/%sModelParams.bin")).exists())) {
            return null;
        }
        //Load parameters from disk:
        INDArray newParams;
        DataInputStream dis = new DataInputStream(new FileInputStream(getPath(modelNamePrefix, "/%sModelParams.bin")));
        newParams = Nd4j.read(dis);

        //Load network configuration from disk:
        MultiLayerConfiguration confFromJson =
                MultiLayerConfiguration.fromJson(FileUtils.readFileToString(new File(getPath(modelNamePrefix, "/%sModelConf.json")), Charset.defaultCharset()));

        //Create a MultiLayerNetwork from the saved configuration and parameters
        net = new MultiLayerNetwork(confFromJson);
        net.init();
        net.setParameters(newParams);

        return net;
    }


    private String getPath(String modelNamePrefix, String format) {
        return modelPath + String.format(format, modelNamePrefix);
    }


}
