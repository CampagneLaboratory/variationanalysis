package org.campagnelab.dl.somatic.mappers;

import org.campagnelab.dl.varanalysis.protobuf.BaseInformationRecords;

import java.util.Arrays;
import java.util.List;
import java.util.Properties;
import java.util.function.Function;

/**
 * Created by rct66 on 12/29/16.
 */
public class DensityMapperCapped extends DensityMapper {


    int linearBinMin;
    int linearBinMax;

    int numLinearBins;

    public DensityMapperCapped(String name1, int linearBinMin, int linearBinMax, Properties sbiProperties,
                               Function<BaseInformationRecords.BaseInformationOrBuilder, List<BaseInformationRecords.NumberWithFrequency>> recordToValues) {
        super(name1, 1, sbiProperties, recordToValues, Integer::floatValue);
        //now cap bins:
        this.numBins = linearBinMax - linearBinMin;
        this.linearBinMax = linearBinMax;
        this.linearBinMin = linearBinMin;
        bins = new float[this.numBins];
        this.recordToValues = recordToValues;
        this.binWidth = 1;
    }

    public DensityMapperCapped(String name1, String name2, int linearBinMin, int linearBinMax, Properties sbiProperties,
                               Function<BaseInformationRecords.BaseInformationOrBuilder, List<BaseInformationRecords.NumberWithFrequency>> recordToValues
    ) {

        super(name1, name2, 1, sbiProperties, recordToValues, Integer::floatValue);
        //now cap bins:
        this.numBins = linearBinMax - linearBinMin;
        this.linearBinMax = linearBinMax;
        this.linearBinMin = linearBinMin;
        bins = new float[this.numBins];
        this.recordToValues = recordToValues;
        this.binWidth = 1;
    }


    public void prepareToNormalize(BaseInformationRecords.BaseInformationOrBuilder record, int indexOfRecord) {
        Arrays.fill(bins, 0);
        List<BaseInformationRecords.NumberWithFrequency> listOfValues = recordToValues.apply(record);
        float numElements = 0;
        for (BaseInformationRecords.NumberWithFrequency n : listOfValues) {
            int featureIndex = (int) ((valueFunction.apply(n.getNumber()) - linearBinMin));
            //handle higher than linearMax case, lower than linearMin case
            if (featureIndex >= (numLinearBins) || featureIndex < 0) {
            } else {
                bins[featureIndex] += n.getFrequency();
                numElements += n.getFrequency();
            }
        }
        // normalize the counts to produce a density:
        if (numElements > 0) {
            for (int featureIndex = 0; featureIndex < numBins; featureIndex++) {
                bins[featureIndex] /= numElements;
            }
        }
    }
}
