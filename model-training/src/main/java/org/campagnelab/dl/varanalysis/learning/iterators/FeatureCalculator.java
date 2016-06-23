package org.campagnelab.dl.varanalysis.learning.iterators;

import org.campagnelab.dl.varanalysis.learning.mappers.FeatureMapper;
import org.campagnelab.dl.varanalysis.learning.mappers.LabelMapper;

/**
 * Maps base information to features for maching learning with neural nets.
 * Created by fac2003 on 5/21/16.
 *
 * @author Fabien Campagne
 */
public interface FeatureCalculator extends FeatureMapper, LabelMapper {


}
