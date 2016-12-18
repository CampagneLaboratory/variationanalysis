package org.campagnelab.dl.framework.mappers;


/**
 * Maps base information to features for machine learning with neural nets.
 * Created by fac2003 on 5/21/16.
 *
 * @author Fabien Campagne
 */
public interface FeatureCalculator<T> extends FeatureMapper<T>, LabelMapper<T> {


}
