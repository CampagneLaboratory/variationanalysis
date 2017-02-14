package org.campagnelab.dl.genotype.mappers;

import org.campagnelab.dl.framework.mappers.AbstractFeatureMapper1D;
import org.campagnelab.dl.varanalysis.protobuf.BaseInformationRecords;

/**
 * Indicate if a count corresponds to a SNP, insertion, deletion, or ref
 */

public class VariantTypeMapper extends AbstractFeatureMapper1D<BaseInformationRecords.BaseInformationOrBuilder> {


    int sampleIndex;
    int genotypeIndex;



    public VariantTypeMapper(int sample, int genotype) {
        this.sampleIndex = sample;
        this.genotypeIndex = genotype;
    }

    @Override
    public int numberOfFeatures() {
        return 4;
    }


    public void prepareToNormalize(BaseInformationRecords.BaseInformationOrBuilder record, int indexOfRecord) {
        //normalization will take place after concating this mapper with other count mappers
    }

    public float produceFeature(BaseInformationRecords.BaseInformationOrBuilder record, int featureIndex) {
        return produceFeatureInternal(record, featureIndex);
    }

    @Override
    public String getFeatureName(int featureIndex) {
        return "variantType" + featureIndex;
    }

    private float normalize(float value, int normalizationFactor) {
        //we are not normalizing here
        return value;
    }


    public float produceFeatureInternal(BaseInformationRecords.BaseInformationOrBuilder record, int featureIndex) {
        assert (featureIndex >= 0 && featureIndex <= 3) : "There are only 4 types of features this mapper corresponds to (SNP, insertion, deletion, ref)";
        int featureProduced = -1;
        if (sampleIndex < record.getSamplesCount()) {
            final BaseInformationRecords.SampleInfo sample = record.getSamples(sampleIndex);
            if (genotypeIndex < sample.getCountsCount()) {
                BaseInformationRecords.CountInfo genoInfo = sample.getCounts(genotypeIndex);
                if (sample.getIsVariant()) {
                    if (genoInfo.getToSequence().length() == 1) {
                        featureProduced = 0;
                    } else {
                        if (genoInfo.getToSequence().contains("-")) {
                            featureProduced = 2;
                        } else {
                            featureProduced = 1;
                        }
                    }
                } else {
                    featureProduced = 3;
                }

            }
        }
        return featureProduced != -1 ? (featureProduced == featureIndex ? 1F : 0F) : -1F;
    }

}
