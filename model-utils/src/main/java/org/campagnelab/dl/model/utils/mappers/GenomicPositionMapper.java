package org.campagnelab.dl.model.utils.mappers;

import org.campagnelab.dl.varanalysis.protobuf.BaseInformationRecords;
import org.nd4j.linalg.api.ndarray.INDArray;

/**
 * Encode a genomic position. Encode both chromosome and position with one hot encoding.
 * Created by fac2003 on 7/12/16.
 */
public class GenomicPositionMapper implements FeatureMapper {
    private static final int NUM_POSITION_BITS = 32;
    private static final int NUM_CHROMOSOMES = 100;
    private final OneHotFeatureMapper chromosomeMapper;
    private final OneHotFeatureMapper positionMapper;
    private ConcatFeatureMapper delegate;

    public GenomicPositionMapper() {
        this.chromosomeMapper = new OneHotFeatureMapper() {

            @Override
            public int getIntegerValue(BaseInformationRecords.BaseInformationOrBuilder record) {
                return record.getReferenceIndex();
            }
        };
        this.positionMapper = new OneHotFeatureMapper() {
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
}
