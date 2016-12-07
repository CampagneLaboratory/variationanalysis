package org.campagnelab.dl.genotype.mappers;

import org.campagnelab.dl.framework.mappers.FeatureMapper;
import org.campagnelab.dl.framework.mappers.FeatureNameMapper;
import org.campagnelab.dl.somatic.genotypes.BaseGenotypeCountFactory;
import org.campagnelab.dl.somatic.genotypes.GenotypeCountFactory;
import org.campagnelab.dl.somatic.mappers.AbstractFeatureMapper;
import org.campagnelab.dl.somatic.mappers.GenotypeCount;
import org.campagnelab.dl.varanalysis.protobuf.BaseInformationRecords;
import org.nd4j.linalg.api.ndarray.INDArray;

/**
 * This is a normalizing mapper. It accepts a single delegate, and normalizes it by dividing all the features by magnitude of the biggest feature.
 *
 * @author Remi Torracinta
 */

public class AbsoluteNormalizeMapper extends AbstractFeatureMapperStripped<BaseInformationRecords.BaseInformationOrBuilder> {

    FeatureNameMapper delegate;
    float absMax=Float.NEGATIVE_INFINITY;

    public AbsoluteNormalizeMapper(FeatureNameMapper delegate) {
        this.delegate = delegate;
    }

    @Override
    public int numberOfFeatures() {
        return delegate.numberOfFeatures();
    }

    public void prepareToNormalize(BaseInformationRecords.BaseInformationOrBuilder record, int indexOfRecord) {
        for (int i = 0; i < numberOfFeatures(); i++) {
            absMax = Math.max(Math.abs(delegate.produceFeature(record, i)), absMax);
        }
    }

    public float produceFeature(BaseInformationRecords.BaseInformationOrBuilder record, int featureIndex) {
        return normalize(produceFeatureInternal(record, featureIndex), absMax);
    }

    @Override
    public String getFeatureName(int featureIndex) {
        return delegate.getFeatureName(featureIndex);
    }


    private float normalize(float value, float normalizationFactor) {
        if (normalizationFactor == 0) {
            return 0;
        }
        float normalized = value / normalizationFactor;
        assert normalized >= -1 && normalized <= 1 : "value must be normalized: " + normalized;
        return normalized;
    }


    private float produceFeatureInternal(BaseInformationRecords.BaseInformationOrBuilder record, int featureIndex) {
        return delegate.produceFeature(record, featureIndex);
    }

}


