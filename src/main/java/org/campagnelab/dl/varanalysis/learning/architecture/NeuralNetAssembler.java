package org.campagnelab.dl.varanalysis.learning.architecture;

import org.deeplearning4j.nn.conf.MultiLayerConfiguration;
import org.deeplearning4j.nn.weights.WeightInit;

/**
 * An interface for classes that assemble neural network configuration according to some architecture.
 * Created by fac2003 on 6/10/16.
 *
 * @author Fabien Campagne
 */
public interface NeuralNetAssembler {
    public MultiLayerConfiguration createNetwork();

    public void setWeightInitialization(WeightInit init);
    public void setLearningRate(double learningRate);

    public void setSeed(int seed);

  //  public void setInputs(int numInputs, int numOutputs, int[] numHiddenNodes);

    public void setNumInputs(int numInputs);

    public void setNumOutputs(int numOutputs);

    public void setNumHiddenNodes(int numHiddenNodes);

    public void setRegularization(boolean regularization);

    public void setRegularizationRate(double regularizationRate);
}
