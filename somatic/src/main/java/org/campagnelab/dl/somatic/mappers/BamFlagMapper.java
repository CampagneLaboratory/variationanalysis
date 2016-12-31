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


}
