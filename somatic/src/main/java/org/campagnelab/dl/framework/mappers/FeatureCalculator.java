package org.campagnelab.dl.framework.mappers;


import org.campagnelab.dl.framework.mappers.FeatureMapper;
import org.campagnelab.dl.framework.mappers.LabelMapper;

/**
 * Maps base information to features for maching learning with neural nets.
 * Created by fac2003 on 5/21/16.
 *
 * @author Fabien Campagne
 */
public interface FeatureCalculator<T> extends FeatureMapper<T>, LabelMapper<T> {


}
