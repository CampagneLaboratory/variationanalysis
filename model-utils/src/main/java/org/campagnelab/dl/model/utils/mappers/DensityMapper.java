package org.campagnelab.dl.model.utils.mappers;

import org.campagnelab.dl.varanalysis.protobuf.BaseInformationRecords;
import org.nd4j.linalg.api.ndarray.INDArray;

import java.util.Arrays;
import java.util.List;
import java.util.Properties;
import java.util.function.Function;

/**
 * Produces feature that represent a density of values for a given number of bins..
 * Created by fac2003 on 10/21/16.
 */
public class DensityMapper implements FeatureMapper, EfficientFeatureMapper, FeatureNameMapper {
    private final Function<BaseInformationRecords.BaseInformationOrBuilder, List<BaseInformationRecords.NumberWithFrequency>> recordToValues;
    private final float minValue;
    private final float maxValue;
    private final float binWidth;
    private final String name;
    int numBins = 10;
    float[] bins;
    private int[] indices;

    public DensityMapper(String name, int numBins, Properties sbiProperties,
                         Function<BaseInformationRecords.BaseInformationOrBuilder, List<BaseInformationRecords.NumberWithFrequency>> recordToValues) {

        if (!propertiesPresent(sbiProperties, "stats." + name)) {
            throw new UnsupportedOperationException("The sbip file does not contain the statistics for " +name+  " (stats."+name+"+.min and stats."+name+".max)");
        }
        this.minValue = getMin(sbiProperties, "stats." + name);
        this.maxValue = getMax(sbiProperties, "stats." + name);


        this.name = name;
        this.numBins = numBins;
        bins = new float[numBins];
        this.recordToValues = recordToValues;
        this.binWidth = (maxValue - minValue) / numBins;
    }


    @Override
    public int numberOfFeatures() {
        return numBins;
    }

    @Override
    public void prepareToNormalize(BaseInformationRecords.BaseInformationOrBuilder record, int indexOfRecord) {
        Arrays.fill(bins, 0);
        List<BaseInformationRecords.NumberWithFrequency> listOfValues = recordToValues.apply(record);
        float numElements = 0;
        for (BaseInformationRecords.NumberWithFrequency n : listOfValues) {
            int featureIndex = (int) ((n.getNumber() - minValue) / binWidth);
            if (featureIndex < 0 || featureIndex >= numBins) {
                //ignore points outside of min-max
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

    @Override
    public void mapFeatures(BaseInformationRecords.BaseInformationOrBuilder record, INDArray inputs, int indexOfRecord) {
        mapFeatures(record, bins, 0, indexOfRecord);
        for (int featureIndex = 0; featureIndex < numberOfFeatures(); featureIndex++) {
            indices[1] = featureIndex;
            inputs.putScalar(indices, bins[featureIndex]);
        }
    }

    @Override
    public float produceFeature(BaseInformationRecords.BaseInformationOrBuilder record, int featureIndex) {
        return bins[featureIndex];
    }

    @Override
    public void mapFeatures(BaseInformationRecords.BaseInformationOrBuilder record, float[] inputs, int offset, int indexOfRecord) {
        // do not copy if inputs is bins (call from mapFeatures above)
        if (inputs != bins) {
            System.arraycopy(bins, 0, inputs, offset, numBins);
        }
    }

    @Override
    public String getFeatureName(int featureIndex) {
        float binMin = 0;
        float binMax = 0;
        for (int i = 0; i < numBins; i++) {
            if (i < featureIndex) {
                binMin += binWidth;
            }
            if (i <= featureIndex) {
                binMax += binWidth;
            }
        }
        return String.format("density_%s_%d_%d", name, binMin, binMax);
    }


    private boolean propertiesPresent(Properties sbiProperties, String s) {
        return sbiProperties.containsKey(s + ".min") || !sbiProperties.containsKey(s + ".max");
    }

    private float getMin(Properties sbiProperties, String propertyName) {
        return Float.parseFloat(sbiProperties.getProperty(propertyName + ".min"));
    }

    private float getMax(Properties sbiProperties, String propertyName) {
        return Float.parseFloat(sbiProperties.getProperty(propertyName + ".max"));
    }

}
