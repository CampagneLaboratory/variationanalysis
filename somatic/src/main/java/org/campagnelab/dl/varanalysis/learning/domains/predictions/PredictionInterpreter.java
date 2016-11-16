package org.campagnelab.dl.varanalysis.learning.domains.predictions;

import org.nd4j.linalg.api.ndarray.INDArray;

/**
 * Instances of this class interpret the numeric predictions of a model into a form understandable
 * by humans.
 * @author  Fabien Campagne
 * Created by fac2003 on 11/12/16.
 */
public interface PredictionInterpreter<RecordType, PredictionType extends Prediction> {
    PredictionType interpret(RecordType record, INDArray output);
}
