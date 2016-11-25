package org.campagnelab.dl.somatic.mappers;

import org.campagnelab.dl.framework.iterators.ConcatFeatureMapper;
import org.campagnelab.dl.framework.mappers.FeatureNameMapper;
import org.campagnelab.dl.framework.mappers.NoMaskFeatureMapper;
import org.campagnelab.dl.framework.mappers.OneHotHashModuloMapper;
import org.campagnelab.dl.varanalysis.protobuf.BaseInformationRecords;
import org.nd4j.linalg.api.ndarray.INDArray;

/**
 * Encode a genomic position. Encode both chromosome and position with one hot encoding hash modulo.
 * Created by fac2003 on 7/12/16.
 */
public class GenomicPositionDiscreteMapper extends NoMaskFeatureMapper<BaseInformationRecords.BaseInformationOrBuilder>
        implements FeatureNameMapper<BaseInformationRecords.BaseInformationOrBuilder> {
    private ConcatFeatureMapper delegate;
    int numChromosomeFeatures = 50;
    int numPositionFeatures = 25000;

    public GenomicPositionDiscreteMapper() {

        OneHotHashModuloMapper<BaseInformationRecords.BaseInformationOrBuilder> chromosomeMapper =
                new OneHotHashModuloMapper<>(numChromosomeFeatures,
                        BaseInformationRecords.BaseInformationOrBuilder::getReferenceId);

        OneHotHashModuloMapper<BaseInformationRecords.BaseInformationOrBuilder> positionMapper =
                new OneHotHashModuloMapper<BaseInformationRecords.BaseInformationOrBuilder>(
                        // we bin position in 10,000 base windows:
                        numPositionFeatures, record -> {
                    return record.getPosition() / 10000;
                }
                );
        this.delegate = new ConcatFeatureMapper(chromosomeMapper, positionMapper);
    }

    @Override
    public int numberOfFeatures() {
        return delegate.numberOfFeatures();

    }

    @Override
    public void prepareToNormalize(BaseInformationRecords.BaseInformationOrBuilder record, int indexOfRecord) {
        delegate.prepareToNormalize(record, indexOfRecord);
    }

    int[] indices = new int[]{0, 0};

    @Override
    public void mapFeatures(BaseInformationRecords.BaseInformationOrBuilder record, INDArray inputs, int indexOfRecord) {
        delegate.mapFeatures(record, inputs, indexOfRecord);
    }

    @Override
    public float produceFeature(BaseInformationRecords.BaseInformationOrBuilder record, int featureIndex) {
        return delegate.produceFeature(record, featureIndex);
    }

    @Override
    public String getFeatureName(int featureIndex) {
        if (featureIndex < numChromosomeFeatures) {
            return "chromosomeHash" + featureIndex;
        } else {
            return "positionWindowHash" + (featureIndex - numChromosomeFeatures);
        }
    }
}
