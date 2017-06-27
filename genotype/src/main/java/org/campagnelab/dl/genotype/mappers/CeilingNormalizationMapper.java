package org.campagnelab.dl.genotype.mappers;

import org.campagnelab.dl.framework.mappers.AbstractFeatureMapper1D;
import org.campagnelab.dl.framework.mappers.FeatureNameMapper;

/**
 * This is a normalizing mapper. It accepts a single delegate, and normalizes the value v of the delegates
 * between 0 a and 1, with 1 the maximum value when v has reached the ceiling. Formally, return min(ceiling,v)/ceiling
 * for positive values and max(ceiling,v)/ceiling for negative values of v.
 *
 * @author Fabien Campagne
 */

public class CeilingNormalizationMapper<RecordType> extends AbstractFeatureMapper1D<RecordType> {

    FeatureNameMapper<RecordType> delegate;
    private float ceiling;

    public CeilingNormalizationMapper(FeatureNameMapper delegate, float ceiling) {
        this.delegate = delegate;
        this.ceiling=ceiling;
    }

    @Override
    public int numberOfFeatures() {
        return delegate.numberOfFeatures();
    }

    public void prepareToNormalize(RecordType record, int indexOfRecord) {
    }

    public float produceFeature(RecordType record, int featureIndex) {
        return normalize(produceFeatureInternal(record, featureIndex));
    }

    @Override
    public String getFeatureName(int featureIndex) {
        return delegate.getFeatureName(featureIndex);
    }


    private float normalize(float value) {
        if (value < 0){
            return Math.max(ceiling, value)/ceiling;
        } else {
            return Math.min(ceiling, value)/ceiling;
        }
    }


    private float produceFeatureInternal(RecordType record, int featureIndex) {
        return delegate.produceFeature(record, featureIndex);
    }

}


