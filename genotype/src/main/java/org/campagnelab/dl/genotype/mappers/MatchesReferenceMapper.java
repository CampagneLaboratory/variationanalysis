package org.campagnelab.dl.genotype.mappers;

import org.campagnelab.dl.framework.mappers.AbstractFeatureMapper1D;
import org.campagnelab.dl.varanalysis.protobuf.BaseInformationRecords;

/**
 * Indicate if a count matches the reference.
 */

public class MatchesReferenceMapper extends AbstractFeatureMapper1D<BaseInformationRecords.BaseInformationOrBuilder> {


    int sampleIndex;
    int genotypeIndex;



    public MatchesReferenceMapper(int sample, int genotype) {
        this.sampleIndex = sample;
        this.genotypeIndex = genotype;
    }

    @Override
    public int numberOfFeatures() {
        return 2;
    }


    public void prepareToNormalize(BaseInformationRecords.BaseInformationOrBuilder record, int indexOfRecord) {
        //normalization will take place after concating this mapper with other count mappers
    }

    public float produceFeature(BaseInformationRecords.BaseInformationOrBuilder record, int featureIndex) {
        return produceFeatureInternal(record, featureIndex);
    }

    @Override
    public String getFeatureName(int featureIndex) {
    return "matchesReference"+featureIndex;
    }

    private float normalize(float value, int normalizationFactor) {
        //we are not normalizing here
        return value;
    }


    public float produceFeatureInternal(BaseInformationRecords.BaseInformationOrBuilder record, int featureIndex) {
        assert (featureIndex >= 0 && featureIndex <= 1) : "This mapper only outputs 1 feature corresponding to one base count";
        BaseInformationRecords.CountInfo genoInfo;
        if (sampleIndex < record.getSamplesCount()) {
            final BaseInformationRecords.SampleInfo sample = record.getSamples(sampleIndex);
            if (genotypeIndex < sample.getCountsCount()) {
                genoInfo = sample.getCounts(genotypeIndex);

                if (featureIndex==0) {
                    return genoInfo.getMatchesReference()?1:0;
                }else{
                    return !genoInfo.getMatchesReference()?1:0;
                }
            }
        }
        return -1;
    }

}
