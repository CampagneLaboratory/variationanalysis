package org.campagnelab.dl.varanalysis.learning.models;

import org.apache.commons.io.FilenameUtils;
import org.deeplearning4j.earlystopping.EarlyStoppingModelSaver;
import org.deeplearning4j.nn.graph.ComputationGraph;
import org.deeplearning4j.util.ModelSerializer;

import java.io.IOException;
import java.nio.charset.Charset;

/**
 * Save a computational graph to disk. Adapted from DL4J model saver, but supports different model prefixes (e.g., best, latest, 1- for epochs,
 * calibrated, etc.)
 *
 * @author Fabien Campagne
 */
public class ComputationGraphSaver implements EarlyStoppingModelSaver<ComputationGraph> {

    private static final String bestFileName = "bestModel.bin";
    private static final String latestFileName = "latestModel.bin";
    private String directory;
    private Charset encoding;

    /**
     * Constructor that uses default character set for configuration (json) encoding
     *
     * @param directory Directory to save networks
     */
    public ComputationGraphSaver(String directory) {
        this(directory, Charset.defaultCharset());
    }

    /**
     * @param directory Directory to save networks
     * @param encoding  Character encoding for configuration (json)
     */
    public ComputationGraphSaver(String directory, Charset encoding) {
        this.directory = directory;
        this.encoding = encoding;
    }


    public void saveModel(ComputationGraph net, String prefix) throws IOException {
        String confOut = FilenameUtils.concat(directory, prefix + "Model.bin");
        save(net, confOut);
    }


    @Override
    public void saveBestModel(ComputationGraph net, double score) throws IOException {
        saveModel(net, "best");
    }

    @Override
    public void saveLatestModel(ComputationGraph net, double score) throws IOException {
        saveModel(net, "latest");

    }

    @Override
    public ComputationGraph getBestModel() throws IOException {
        String confOut = FilenameUtils.concat(directory, bestFileName);
        return load(confOut);
    }

    @Override
    public ComputationGraph getLatestModel() throws IOException {
        String confOut = FilenameUtils.concat(directory, latestFileName);
        return load(confOut);
    }

    private void save(ComputationGraph net, String modelName) throws IOException {
        ModelSerializer.writeModel(net, modelName, true);
    }

    private ComputationGraph load(String modelName) throws IOException {
        return ModelSerializer.restoreComputationGraph(modelName);
    }

    @Override
    public String toString() {
        return "ComputationGraphSaver(dir=" + directory + ")";
    }
}
