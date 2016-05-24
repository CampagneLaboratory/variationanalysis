package org.campagnelab.dl.varanalysis.learning.iterators;

import org.campagnelab.dl.varanalysis.format.PosRecord;
import org.nd4j.linalg.api.ndarray.INDArray;
import org.nd4j.linalg.dataset.DataSet;

/**
 * Maps base information to features for maching learning with neural nets.
 * Created by fac2003 on 5/21/16.
 *
 * @author Fabien Campagne
 */
interface FeatureCalculator extends FeatureMapper, LabelMapper {


}
