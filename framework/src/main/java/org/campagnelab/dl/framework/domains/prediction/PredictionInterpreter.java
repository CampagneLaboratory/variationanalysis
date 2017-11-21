package org.campagnelab.dl.framework.domains.prediction;

import org.nd4j.linalg.api.ndarray.INDArray;

/**
 * Instances of this class interpret the numeric predictions of a model into a form understandable
 * by humans.
 * @author  Fabien Campagne
 * Created by fac2003 on 11/12/16.
 */
public interface PredictionInterpreter<RecordType, PredictionType extends Prediction> {
    /**
     * Interpret a prediction given true labels given by a LabelMapper.
     * @param trueLabels True labels for the model output to be interpreted.
     * @param output   Model output to be interpreted.
     * @param exampleIndex Index of the example being predicted, in a mini-batch.
     * @return Interpreted prediction.
     */
    PredictionType interpret(INDArray trueLabels, INDArray output, int exampleIndex);

    /**
     * Interpret a prediction given a record and model outputs.
     * @param record The record, which can be mapped to true labels with a labelMapper.
     * @param output Model outputs. In the order defined by the computation graph that produced them.
     * @return Interpreted prediction.
     */
    PredictionType interpret(RecordType record, INDArray output);
}
