package org.campagnelab.dl.framework.models;

import org.apache.commons.io.FilenameUtils;
import org.deeplearning4j.earlystopping.EarlyStoppingModelSaver;
import org.deeplearning4j.nn.multilayer.MultiLayerNetwork;
import org.deeplearning4j.util.ModelSerializer;

import java.io.IOException;
import java.nio.charset.Charset;

/** Save a model to disk. Adapted from DL4J model saver, but supports different model prefixes (e.g., best, latest, 1- for epochs,
 * calibrated, etc.)
 * @author Fabien Campagne
 */
public class ModelSaver implements EarlyStoppingModelSaver<MultiLayerNetwork> {

    private static final String bestFileName = "bestModel.bin";
    private static final String latestFileName = "latestModel.bin";
    private String directory;
    private Charset encoding;

    /**Constructor that uses default character set for configuration (json) encoding
     * @param directory Directory to save networks
     */
    public ModelSaver(String directory) {
        this(directory, Charset.defaultCharset());
    }

    /**
     * @param directory Directory to save networks
     * @param encoding Character encoding for configuration (json)
     */
    public ModelSaver(String directory, Charset encoding){
        this.directory = directory;
        this.encoding = encoding;
    }

    @Override
    public void saveBestModel(MultiLayerNetwork net, double score) throws IOException {
        saveModel(net,"best",score);
    }

    public void saveModel(MultiLayerNetwork net,String prefix, double score) throws IOException {
        String confOut = FilenameUtils.concat(directory,prefix+"Model.bin");
        save(net,confOut);
    }

    @Override
    public void saveLatestModel(MultiLayerNetwork net, double score) throws IOException {
      saveModel(net,"latest",score);
    }

    @Override
    public MultiLayerNetwork getBestModel() throws IOException {
        String confOut = FilenameUtils.concat(directory, bestFileName);
        return load(confOut);
    }

    @Override
    public MultiLayerNetwork getLatestModel() throws IOException {
        String confOut = FilenameUtils.concat(directory, latestFileName);
        return load(confOut);
    }

    private void save(MultiLayerNetwork net, String modelName) throws IOException{
        ModelSerializer.writeModel(net, modelName, true);
    }

    private MultiLayerNetwork load(String modelName) throws IOException {
        MultiLayerNetwork net = ModelSerializer.restoreMultiLayerNetwork(modelName);
        return net;
    }

    @Override
    public String toString(){
        return "ModelSaver(dir=" + directory + ")";
    }
}
