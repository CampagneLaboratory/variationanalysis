package org.campagnelab.dl.genotype.mappers;

import org.campagnelab.dl.framework.mappers.ConfigurableFeatureMapper;
import org.campagnelab.dl.somatic.mappers.*;
import org.campagnelab.dl.somatic.mappers.functional.TraversalHelper;
import org.campagnelab.dl.varanalysis.protobuf.BaseInformationRecords;
import org.nd4j.linalg.api.ndarray.INDArray;

import java.util.Properties;


/**
 * A somatic feature mapper that reuses GenotypeMapperV28 and adds somatic specific features.
 */
public class SomaticFeatureMapperV26 extends NamingConcatFeatureMapper<BaseInformationRecords.BaseInformationOrBuilder>
        implements ConfigurableFeatureMapper {
    private NamingConcatFeatureMapper delegate;

    private String recordTo(final int contextLength, BaseInformationRecords.BaseInformationOrBuilder record, int countIndex) {
        return MappingFunctions.recordTo(contextLength, record, countIndex);
    }

    /**
     * Configure the feature mapper for a specific set of sbi files. This method accesses the properties of the reader.
     *
     * @param sbiProperties properties from an sbi reader.
     */
    public void configure(Properties sbiProperties) {


        final GenotypeMapperV28 a = new GenotypeMapperV28(0);
        final GenotypeMapperV28 b = new GenotypeMapperV28(1);
        a.configure(sbiProperties);
        b.configure(sbiProperties);
        delegate = new NamingConcatFeatureMapper<BaseInformationRecords.BaseInformationOrBuilder>(
                a,
                b,
                new FractionDifferences4(),
                new MagnitudeFeatures2(),
                new SimpleFeatureCalculator(true)
        );


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

}
