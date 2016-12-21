package org.campagnelab.dl.genotype.predictions;

import org.campagnelab.dl.framework.models.ModelLoader;
import org.campagnelab.goby.algorithmic.dsv.DiscoverVariantPositionData;
import org.campagnelab.goby.algorithmic.dsv.SampleCountInfo;
import org.campagnelab.goby.predictions.GenotypePredictor;
import org.campagnelab.goby.predictions.Predictor;
import org.campagnelab.goby.reads.RandomAccessSequenceInterface;

import java.io.IOException;
import java.util.Properties;

/**
 * This class implements the genotype prediction model expected by Goby 3.2+.  Make sure to move
 * the service declaration (currently in  goby-spi/META-INF.services).
 *
 * @author Fabien Campagne
 *         Created by fac2003 on 11/14/16.
 */
public class DLGenotypePredictor implements GenotypePredictor, Predictor {
    private GenotypeModel model;
    private GenotypePrediction prediction;



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
        model = new GenotypeModel(modelPath, modelPrefix);
    }

    @Override
    public void predict(RandomAccessSequenceInterface genome, String referenceId, SampleCountInfo[] sampleCounts,
                        int referenceIndex, int pos, DiscoverVariantPositionData list, int[] readerIdxs) {
        prediction = model.predictGenotype(genome,
                referenceId,
                sampleCounts,
                referenceIndex, pos,
                list,
                readerIdxs);
    }


    @Override
    public boolean modelIsLoaded() {
        return model != null;
    }

    @Override
    public double probabilityGenotypeIsCalled(int genotypeIndex) {
        throw new NoSuchMethodError();
    }

    @Override
    public double probabilityGenotypeIsNotCalled(int genotypeIndex) {
        throw new NoSuchMethodError();
    }

    @Override
    public String getCalledGenotype() {
        return prediction.predictedGenotype;
    }

    /**
     * The probability of the called genotype, as estimaed by the model.
     *
     * @return a number between 0 and 1 inclusive.
     */
    public double getProbabilityOfCalledGenotype() {
        return prediction.overallProbability;
    }

    @Override
    public boolean trainedForIndels() {
        return false;
    }

    @Override
    public Properties getModelProperties() {

        assert modelIsLoaded() : "You must load a model before you can access model properties.";
        return model.getProperties();
    }

}
