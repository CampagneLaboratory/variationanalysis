package org.campagnelab.dl.genotype.mappers;

import org.campagnelab.dl.framework.mappers.AbstractFeatureMapper1D;
import org.campagnelab.dl.varanalysis.protobuf.BaseInformationRecords;

/**
 * This is a simple feature mapper.
 * Each sampleIndex has the following information:
 * <pre> 69033640	11	false
 * position=14521   referenceIndex=0       isMutated=false
 * sampleIndex 0 counts=[10, 0, 0, 63, 0, 4, 0, 1, 62, 0]
 *
 * Using these data, we can normalize counts and map them
 *
 * Created by rct66 on 11/16/16.
 *
 * @author Remi Torracinta
 */

public class SingleGenoTypeCountMapper extends AbstractFeatureMapper1D<BaseInformationRecords.BaseInformationOrBuilder> {


    int sampleIndex;
    int genotypeIndex;
    boolean getForwardStrand;


    public SingleGenoTypeCountMapper(int sample, int genotype, boolean getForwardStrand) {
        this.sampleIndex = sample;
        this.genotypeIndex = genotype;
        this.getForwardStrand = getForwardStrand;
    }

    @Override
    public int numberOfFeatures() {
        return 1;
    }


    public void prepareToNormalize(BaseInformationRecords.BaseInformationOrBuilder record, int indexOfRecord) {
        //normalization will take place after concating this mapper with other count mappers
    }

    public float produceFeature(BaseInformationRecords.BaseInformationOrBuilder record, int featureIndex) {
        return produceFeatureInternal(record, featureIndex);
    }

    @Override
    public String getFeatureName(int featureIndex) {
        String sampleString = Integer.toString(sampleIndex) + "Sample";
        String genotypeString;
        switch (genotypeIndex) {
            case 0:
                genotypeString = "A";
                break;
            case 1:
                genotypeString = "T";
                break;
            case 2:
                genotypeString = "C";
                break;
            case 3:
                genotypeString = "G";
                break;
            case 4:
                genotypeString = "N";
                break;
            case 5:
                genotypeString = "Indel1";
                break;
            case 6:
                genotypeString = "Indel2";
                break;
            case 7:
                genotypeString = "Indel3";
                break;
            case 8:
                genotypeString = "Indel4";
                break;
            default:
                genotypeString = "Indel5";
        }
        return ("sampleIndex" + sampleString + "Genotype" + genotypeString + (getForwardStrand ? "Forward" : "Backward"));
    }

    private float normalize(float value, int normalizationFactor) {
        //we are not normalizing here
        return value;
    }


    public float produceFeatureInternal(BaseInformationRecords.BaseInformationOrBuilder record, int featureIndex) {
        assert (featureIndex >= 0 && featureIndex < 1) : "This mapper only outputs 1 feature corresponding to one base count";
        BaseInformationRecords.CountInfo genoInfo;
        if (sampleIndex < record.getSamplesCount()) {
            final BaseInformationRecords.SampleInfo sample = record.getSamples(sampleIndex);
            if (genotypeIndex < sample.getCountsCount()) {
                genoInfo = sample.getCounts(genotypeIndex);
                float val = getForwardStrand ? genoInfo.getGenotypeCountForwardStrand() : genoInfo.getGenotypeCountReverseStrand();
                return val;
            }
        }
        return -1;
    }

}
