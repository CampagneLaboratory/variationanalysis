package org.campagnelab.dl.somatic.predictions;

import org.campagnelab.dl.framework.models.ModelLoader;
import org.campagnelab.dl.somatic.utils.ProtoPredictor;
import org.campagnelab.goby.algorithmic.dsv.DiscoverVariantPositionData;
import org.campagnelab.goby.algorithmic.dsv.SampleCountInfo;
import org.campagnelab.goby.predictions.SomaticPredictor;
import org.campagnelab.goby.reads.RandomAccessSequenceInterface;

import java.io.IOException;

/**
 * This class will be moved to the variation project to remove the dependency on model-utils, then loaded
 * with the Java service infrastructure from the model-util jar added in the goby wrapper. Make sure to move
 * the service declaration (currently in  goby-spi/META-INF.services).
 *
 * @author Fabien Campagne
 *         Created by fac2003 on 11/14/16.
 */
public class DLSomaticPredictor implements SomaticPredictor {
    private SomaticModel model;
    private ProtoPredictor.Prediction prediction;

    @Override
    public String getModelPath(String fullMPath) {
        return ModelLoader.getModelPath(fullMPath);
    }

    @Override
    public String getModelPrefix(String fullMPath) {
        return ModelLoader.getModelLabel(fullMPath);
    }

    @Override
    public void loadModel(String modelPath, String modelPrefix) throws IOException {
        model = new SomaticModel(modelPath, modelPrefix);
    }

    @Override
    public void predict(RandomAccessSequenceInterface genome, String referenceId, SampleCountInfo[] sampleCounts,
                        int referenceIndex, int pos, DiscoverVariantPositionData list, int[] readerIdxs) {
        prediction = model.mutPrediction(genome,
                referenceId,
                sampleCounts,
                referenceIndex, pos,
                list,
                readerIdxs);
    }

    @Override
    public double probabilityIsMutated() {

        assert prediction != null;
        return prediction.posProb;
    }

    @Override
    public double probabilityIsNotMutated() {

        assert prediction != null;
        return prediction.negProb;
    }

    @Override
    public float getSomaticFrequency() {

        assert prediction != null;
        return prediction.predictedSomaticFrequency;
    }

    @Override
    public boolean modelIsLoaded() {
        return model != null;
    }

    @Override
    public boolean hasSomaticFrequency() {
        return prediction.hasSomaticFrequency;
    }
}
