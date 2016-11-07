package org.campagnelab.dl.model.utils.mappers;


/**
 * Maps base information to features for maching learning with neural nets.
 * Created by fac2003 on 5/21/16.
 *
 * @author Fabien Campagne
 */
public interface FeatureCalculator<T> extends FeatureMapper<T>, LabelMapper<T> {


}
