package org.campagnelab.dl.genotype.mappers;

import org.campagnelab.dl.framework.mappers.AbstractFeatureMapper1D;
import org.campagnelab.dl.framework.mappers.FeatureNameMapper;
import org.campagnelab.dl.varanalysis.protobuf.BaseInformationRecords;

/**
 * This is a normalizing mapper. It accepts a single delegate, and normalizes it by inversing the values
 *
 * @author Remi Torracinta
 */

public class InverseNormalizationMapper<RecordType> extends AbstractFeatureMapper1D<RecordType> {

    FeatureNameMapper<RecordType> delegate;

    public InverseNormalizationMapper(FeatureNameMapper delegate) {
        this.delegate = delegate;
    }

    @Override
    public int numberOfFeatures() {
        return delegate.numberOfFeatures();
    }

    public void prepareToNormalize(RecordType record, int indexOfRecord) {
    }

    public float produceFeature(RecordType record, int featureIndex) {
        return normalize(produceFeatureInternal(record, featureIndex), 0);
    }

    @Override
    public String getFeatureName(int featureIndex) {
        return delegate.getFeatureName(featureIndex);
    }


    private float normalize(float value, float normalizationFactor) {
        if (value < 0){
            value--;
        } else {
            value++;
        }
        float normalized = 1 / (value);
        assert normalized >= -1 && normalized <= 1 : "value must be normalized: " + normalized;
        return normalized;
    }


    private float produceFeatureInternal(RecordType record, int featureIndex) {
        return delegate.produceFeature(record, featureIndex);
    }

}


