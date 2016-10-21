package org.campagnelab.dl.model.utils.mappers;

import it.unimi.dsi.fastutil.ints.IntArrayList;
import it.unimi.dsi.fastutil.objects.ObjectArrayList;
import org.campagnelab.dl.model.utils.ConfigurableFeatureMapper;
import org.campagnelab.dl.model.utils.mappers.functional.TraversalHelper;
import org.campagnelab.dl.varanalysis.protobuf.BaseInformationRecords;
import org.campagnelab.goby.baseinfo.SequenceBaseInformationReader;
import org.campagnelab.goby.baseinfo.StatAccumulatorNumVariationsInRead;
import org.datavec.api.records.reader.RecordReader;
import org.nd4j.linalg.api.ndarray.INDArray;

import java.util.List;
import java.util.Properties;

/**
 * Same as V18, but adds density features for numVariationsInRead. Starting to use Java8 lambdas to customize generic feature mappers.
 */
public class FeatureMapperV19 extends NamingConcatFeatureMapper implements ConfigurableFeatureMapper {
    NamingConcatFeatureMapper delegate;

    /**
     * Configure the feature mapper for a specific set of sbi files. This method accesses the properties of the reader.
     *
     * @param reader
     */
    public void configure(SequenceBaseInformationReader reader) {
        Properties sbiProperties = reader.getProperties();
        if (!sbiProperties.containsKey("stats.numVariationsInRead.min") || !sbiProperties.containsKey("stats.numVariationsInRead.max")) {
            throw new UnsupportedOperationException("The sbip file does not contain the statistics for numVariationsInRead (stats.numVariationsInRead.min and stats.numVariationsInRead.max)");
        }
        float minNumVariationsInRead = Float.parseFloat(sbiProperties.getProperty("stats.numVariationsInRead.min"));
        float maxNumVariationsInRead = Float.parseFloat(sbiProperties.getProperty("stats.numVariationsInRead.max"));

        delegate = new NamingConcatFeatureMapper(new SimpleFeatureCalculator(true),
                new IndelFeatures(),
                new ReadIndexFeaturesFix(),
                new FractionDifferences4(),
                new MagnitudeFeatures2(),
                new DensityMapper("numVariationsInRead", 10, minNumVariationsInRead, maxNumVariationsInRead, baseInformationOrBuilder ->
                        TraversalHelper.forAllSampleCounts(baseInformationOrBuilder, BaseInformationRecords.CountInfo::getNumVariationsInReadsList)));
       // new DensityMapper("mappingQuality", 25, minNumVariationsInRead, maxNumVariationsInRead, baseInformationOrBuilder ->
        //        TraversalHelper.forAllSampleCounts(baseInformationOrBuilder, BaseInformationRecords.CountInfo::getMap)));

    }


    @Override
    public String getFeatureName(int i) {
        return delegate.getFeatureName(i);
    }

    @Override
    public int numberOfFeatures() {
        return delegate.numberOfFeatures();
    }

    @Override
    public void prepareToNormalize(BaseInformationRecords.BaseInformationOrBuilder record, int indexOfRecord) {
        delegate.prepareToNormalize(record, indexOfRecord);
    }

    @Override
    public void mapFeatures(BaseInformationRecords.BaseInformationOrBuilder record, INDArray inputs, int indexOfRecord) {
        delegate.mapFeatures(record, inputs, indexOfRecord);
    }

    @Override
    public float produceFeature(BaseInformationRecords.BaseInformationOrBuilder record, int featureIndex) {
        return delegate.produceFeature(record, featureIndex);
    }

    @Override
    public void mapFeatures(BaseInformationRecords.BaseInformationOrBuilder record, float[] inputs, int offset, int indexOfRecord) {
        delegate.mapFeatures(record, inputs, offset, indexOfRecord);
    }
}
