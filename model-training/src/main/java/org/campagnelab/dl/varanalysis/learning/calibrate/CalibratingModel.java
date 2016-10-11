package org.campagnelab.dl.varanalysis.learning.calibrate;

import it.unimi.dsi.fastutil.floats.FloatArrayList;
import org.campagnelab.dl.model.utils.mappers.FeatureMapper;
import org.deeplearning4j.nn.api.Layer;
import org.deeplearning4j.nn.multilayer.MultiLayerNetwork;
import org.nd4j.linalg.api.ndarray.INDArray;
import org.nd4j.linalg.factory.Nd4j;

/**
 * Created by fac2003 on 7/17/16.
 */
public class CalibratingModel {
    private final int numPredictiveInputFeatures;
    private final int numCalibratingInputFeatures;
    private final FeatureMapper featureMapper;
    MultiLayerNetwork predictiveModel;
    MultiLayerNetwork calibratingModel;
    private int[] indices = {0, 0};

    public CalibratingModel(MultiLayerNetwork predictiveModel, FeatureMapper featureMapper,
                            MultiLayerNetwork calibratingModel) {
        this.predictiveModel = predictiveModel;
        this.featureMapper = featureMapper;
        this.calibratingModel = calibratingModel;
        this.numPredictiveInputFeatures = featureMapper.numberOfFeatures();
        this.numCalibratingInputFeatures = getModelActivationNumber(predictiveModel, featureMapper);
    }

    public float estimateCalibratedP(INDArray testFeatures) {
        assert testFeatures.data().length() == featureMapper.numberOfFeatures() : "number of features does not match predictive model input length.";
        INDArray inputs = Nd4j.zeros(1, numCalibratingInputFeatures);
        //  INDArray labels = Nd4j.zeros(1, 1);
        // DataSet dataset = new DataSet(inputs,labels);
        int indexOfNewRecordInMinibatch = 0;
        // calculate the model output and use as features for the calibration model:
        FloatArrayList floats = getModelInternalActivations(testFeatures);
        assert floats.size() == numCalibratingInputFeatures : "number of features does not match calibrating model input length.";

        indices[0] = indexOfNewRecordInMinibatch;
        for (int i = 0; i < numCalibratingInputFeatures; i++) {
            indices[1] = i;
            //     dataset.getFeatures().putScalar(indices, floats.getFloat(i));
            inputs.putScalar(indices, floats.getFloat(i));
        }
        calibratingModel.init();
       // System.out.println("size:" + inputs);
        float[] predicted = calibratingModel.output(inputs, false).getRow(0).data().asFloat();
        return predicted[0];
    }

    private FloatArrayList getModelInternalActivations(INDArray testFeatures) {


        FloatArrayList floats = new FloatArrayList();
        predictiveModel.feedForward(testFeatures).stream().forEach(indArray -> floats.addAll(FloatArrayList.wrap(indArray.data().asFloat())));
        return floats;
    }

    private int getModelActivationNumber(MultiLayerNetwork model, FeatureMapper modelFeatureMapper) {
        int numActivations = 0;
        Layer[] layers = model.getLayers();
        int totalNumParams = 0;

        INDArray inputFeatures = Nd4j.zeros(1, modelFeatureMapper.numberOfFeatures());

        int sum = model.feedForward(inputFeatures, false).stream().mapToInt(indArray ->
                indArray.data().asFloat().length).sum();
        System.out.println("Number of activations: " + sum);
        return sum;
    }

}


