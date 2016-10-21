package org.campagnelab.dl.model.utils.mappers;

import it.unimi.dsi.fastutil.objects.ObjectArrayList;
import org.campagnelab.dl.varanalysis.protobuf.BaseInformationRecords;
import org.campagnelab.goby.baseinfo.SequenceBaseInformationReader;
import org.nd4j.linalg.api.ndarray.INDArray;

import java.util.Arrays;
import java.util.List;
import java.util.function.Function;

/**
 * Produces feature that represent a density of values for a given number of bins..
 * Created by fac2003 on 10/21/16.
 */
public class DensityMapper implements FeatureMapper, EfficientFeatureMapper, FeatureNameMapper {
    private final Function<BaseInformationRecords.BaseInformationOrBuilder, List<BaseInformationRecords.NumberWithFrequency>> recordToValues;
    private final float minValue;
    private final float binWidh;
    private final String name;
    int numBins = 10;
    float[] bins;
    private int[] indices;

    public DensityMapper(String name, int numBins, float minValue, float maxValue,
                         Function<BaseInformationRecords.BaseInformationOrBuilder, List<BaseInformationRecords.NumberWithFrequency>> recordToValues) {
        this.name = name;
        this.numBins = numBins;
        bins = new float[numBins];
        this.recordToValues = recordToValues;
        this.minValue = minValue;
        this.binWidh = (maxValue - minValue) / numBins;
    }

    @Override
    public void configure(SequenceBaseInformationReader reader) {
        // configuration is done once and for all before calling the constructor to set minValue and maxValue.
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
            int featureIndex = (int) ((n.getNumber() - minValue) / binWidh);
            if (featureIndex < 0 || featureIndex >= numBins) {
                //ignore points outside of min-max
            } else {
                bins[featureIndex] += n.getFrequency();
                numElements += n.getFrequency();
            }
        }
        // normalize the counts to produce a density:
        for (int featureIndex = 0; featureIndex < numBins; featureIndex++) {
            bins[featureIndex] /= numElements;
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
                binMin += binWidh;
            }
            if (i <= featureIndex) {
                binMax += binWidh;
            }
        }
        return String.format("density_%s_%d_%d", name, binMin, binMax);
    }

    /**
     * Define a Function to reduce a record to a list of NumberWithFrequency found across all samples and counts of these samples.
     * @param baseInformationOrBuilder
     * @param function
     * @return
     */
    public static List<BaseInformationRecords.NumberWithFrequency> forAllCounts(BaseInformationRecords.BaseInformationOrBuilder baseInformationOrBuilder,
                                                                                Function<BaseInformationRecords.CountInfo,List<BaseInformationRecords.NumberWithFrequency>> function) {
        List<BaseInformationRecords.NumberWithFrequency> list = new ObjectArrayList<>();

        baseInformationOrBuilder.getSamplesList().forEach(
                sampleInfo -> {
                    sampleInfo.getCountsList().forEach(
                          countInfo -> list.addAll(function.apply(countInfo))
                    );
                }
        );
        return list;
    }
}
