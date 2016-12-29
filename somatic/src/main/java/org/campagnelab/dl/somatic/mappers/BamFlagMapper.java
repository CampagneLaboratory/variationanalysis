package org.campagnelab.dl.somatic.mappers;

import org.campagnelab.dl.framework.mappers.AbstractFeatureMapper1D;
import org.campagnelab.dl.varanalysis.protobuf.BaseInformationRecords;

import java.util.List;

/**
 * Produces feature that represent a density of values for a given number of bins..
 * Created by fac2003 on 10/21/16.
 */
public class BamFlagMapper extends AbstractFeatureMapper1D<BaseInformationRecords.BaseInformationOrBuilder> {

    final static int BOOL_NUM = 12;

    int genotypeIndex;
    int alignmentCount = 0;
    int[] propCounts;
    double[] propFractions;
    int sampleIndex;

    private BamFlagMapper() {
    }

    public BamFlagMapper(int sample, int genotypeIndex) {
        this.sampleIndex = sample;
        this.genotypeIndex = genotypeIndex;
    }

    @Override
    public int numberOfFeatures() {
        return BOOL_NUM;
    }

    @Override
    public void prepareToNormalize(BaseInformationRecords.BaseInformationOrBuilder record, int indexOfRecord) {
        propCounts = new int[BOOL_NUM];
        propFractions = new double[BOOL_NUM];
        List<BaseInformationRecords.NumberWithFrequency> bamFlagFreqs = record.getSamples(sampleIndex).getCounts(genotypeIndex).getPairFlagsList();
        for (BaseInformationRecords.NumberWithFrequency bamFlagFreq : bamFlagFreqs) {
            boolean[] decoded = decodeProps(bamFlagFreq.getNumber());
            for (int i = 0; i < BOOL_NUM; i++) {
                if (decoded[i]) {
                    propCounts[i] += bamFlagFreq.getFrequency();
                }
            }
            alignmentCount += bamFlagFreq.getFrequency();
        }
        for (int i = 0; i < BOOL_NUM; i++) {
            if (alignmentCount != 0) {
                propFractions[i] = propCounts[i] / (double) (alignmentCount);
            } else {
                propFractions[i] = 0;
            }
        }

    }


    boolean[] decodeProps(int bamFlag) {
        boolean[] decoded = new boolean[BOOL_NUM];
        for (int i = 0; i < BOOL_NUM; i++) {
            if (getBit(bamFlag, i) == 1) {
                decoded[i] = true;
            }
        }
        return decoded;
    }

    @Override
    public float produceFeature(BaseInformationRecords.BaseInformationOrBuilder record, int featureIndex) {
        return (float) propFractions[featureIndex];
    }


    @Override
    public String getFeatureName(int featureIndex) {
        switch (featureIndex) {
            case 0:
                return "template having multiple segments in sequencing";
            case 1:
                return "each segment properly aligned according to the aligner";
            case 2:
                return "segment unmapped";
            case 3:
                return "next segment in the template unmapped";
            case 4:
                return "SEQ being reverse complemented";
            case 5:
                return "SEQ of the next segment in the template being reversed";
            case 6:
                return "the first segment in the template";
            case 7:
                return "the last segment in the template";
            case 8:
                return "secondary alignment";
            case 9:
                return "not passing quality controls";
            case 10:
                return "PCR or optical duplicate";
            case 11:
                return "supplementary alignment";
            default:
                return "unknown flag";
        }
    }

    static int getBit(int n, int k) {
        return (n >> k) & 1;
    }

}
