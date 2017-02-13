package org.campagnelab.dl.somatic.mappers;

import org.campagnelab.dl.framework.mappers.AbstractFeatureMapper1D;
import org.campagnelab.dl.varanalysis.protobuf.BaseInformationRecords;

import java.util.Arrays;
import java.util.List;
import java.util.function.Function;

/**
 * Map binary flags to a specific number of bins, as a histogram (each bin indicates the frequency of finding the bit=1).
 */
public class BamFlagMapper extends AbstractFeatureMapper1D<BaseInformationRecords.BaseInformationOrBuilder> {


    private int BOOL_NUM = 12;
    private final Function<BaseInformationRecords.BaseInformationOrBuilder, List<BaseInformationRecords.NumberWithFrequency>> function;
    private final int[] propCounts;
    double[] propFractions;
    /**
     * * @param sampleIndex index of the sample to get the list from.
     *
     * @param genotypeIndex index of the genotype to get the list from.
     */
    public BamFlagMapper(final int sampleIndex, final int genotypeIndex) {
        this(12, baseInformationOrBuilder -> {
            return baseInformationOrBuilder.getSamples(sampleIndex).getCounts(genotypeIndex).getPairFlagsList();
        });
    }

    /**
     * Construct a mapper for lists of boolean flag integers stored in NumberWithFrequency.
     *
     * @param flagSize Number of bits to map in the flag (lower bits are extracted).
     * @param function function that extracts a list of number from a record.
     */
    public BamFlagMapper(int flagSize,
                         Function<BaseInformationRecords.BaseInformationOrBuilder,
                                 List<BaseInformationRecords.NumberWithFrequency>> function) {
        this.BOOL_NUM = flagSize;
        this.function = function;
        propCounts = new int[BOOL_NUM];
        propFractions = new double[BOOL_NUM];
    }

    @Override
    public int numberOfFeatures() {
        return BOOL_NUM;
    }

    @Override
    public void prepareToNormalize(BaseInformationRecords.BaseInformationOrBuilder record, int indexOfRecord) {
        Arrays.fill(propCounts, 0);
        Arrays.fill(propFractions, 0);

        List<BaseInformationRecords.NumberWithFrequency> bamFlagFreqs = function.apply(record);
        int overallSum = 0;
        for (BaseInformationRecords.NumberWithFrequency bamFlagFreq : bamFlagFreqs) {
            int flags = bamFlagFreq.getNumber();
            int frequency = bamFlagFreq.getFrequency();
            // reverser order here so bits and histogram are in the same order:
            for (int i = BOOL_NUM; i >= 0; i--) {
                if (getBit(flags, i) == 1) {
                    propCounts[i] += frequency;
                    overallSum += frequency;
                }
            }
            //  alignmentCount += frequency;
        }
        for (int i = 0; i < BOOL_NUM; i++) {
            if (overallSum != 0) {
                propFractions[i] = propCounts[i] / (double) (overallSum);
            } else {
                propFractions[i] = 0;
            }
        }

    }

    private int getBit(int n, int k) {
        return (n >> k) & 1;
    }

    @Override
    public float produceFeature(BaseInformationRecords.BaseInformationOrBuilder record, int featureIndex) {
        return (float) propFractions[featureIndex];
    }


    @Override
    public String getFeatureName(int featureIndex) {
        switch (featureIndex) {
            case 0:
                return "p_multiple_segs";
            case 1:
                return "p_allsegs_aligned";
            case 2:
                return "p_seg_unmapped";
            case 3:
                return "p_nextseg_unmapped";
            case 4:
                return "p_SEQ_rev_comp";
            case 5:
                return "p_nextSEQ_rev_comp";
            case 6:
                return "p_first_seg";
            case 7:
                return "p_last_seg";
            case 8:
                return "p_secondary_align";
            case 9:
                return "p_fail_qualControl";
            case 10:
                return "p_PCR||Opt_duplicate";
            case 11:
                return "p_supplement_align";
            default:
                return "p_unknown_flag";
        }
    }


}
