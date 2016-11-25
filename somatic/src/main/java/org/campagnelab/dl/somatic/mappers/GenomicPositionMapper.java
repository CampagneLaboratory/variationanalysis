package org.campagnelab.dl.somatic.mappers;

import org.campagnelab.dl.framework.iterators.ConcatFeatureMapper;
import org.campagnelab.dl.framework.mappers.FeatureNameMapper;
import org.campagnelab.dl.framework.mappers.NoMaskFeatureMapper;
import org.campagnelab.dl.varanalysis.protobuf.BaseInformationRecords;
import org.nd4j.linalg.api.ndarray.INDArray;

/**
 * Encode a genomic position. Encode both chromosome and position with one hot encoding.
 * Created by fac2003 on 7/12/16.
 */
public class GenomicPositionMapper extends NoMaskFeatureMapper<BaseInformationRecords.BaseInformationOrBuilder>
implements FeatureNameMapper<BaseInformationRecords.BaseInformationOrBuilder>{
    private static final int NUM_POSITION_BITS = 32;
    private static final int NUM_CHROMOSOMES = 100;
    private final BinaryFeatureMapper chromosomeMapper;
    private final BinaryFeatureMapper positionMapper;
    private ConcatFeatureMapper delegate;

    public GenomicPositionMapper() {
        this.chromosomeMapper = new BinaryFeatureMapper() {

            @Override
            public int getIntegerValue(BaseInformationRecords.BaseInformationOrBuilder record) {
                return record.getReferenceIndex();
            }
        };
        this.positionMapper = new BinaryFeatureMapper() {
            @Override
            public int getIntegerValue(BaseInformationRecords.BaseInformationOrBuilder record) {
                return record.getPosition();
            }
        };
        this.delegate = new ConcatFeatureMapper(chromosomeMapper,positionMapper);
    }

    @Override
    public int numberOfFeatures() {
        return chromosomeMapper.numberOfFeatures() + positionMapper.numberOfFeatures();

    }

    @Override
    public void prepareToNormalize(BaseInformationRecords.BaseInformationOrBuilder record, int indexOfRecord) {

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
        return "genomicPosition"+featureIndex;
    }
}
