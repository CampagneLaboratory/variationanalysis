/*
package org.campagnelab.dl.model.utils.mappers;

import it.unimi.dsi.fastutil.objects.ObjectArrayList;
import org.campagnelab.dl.varanalysis.protobuf.BaseInformationRecords;
import org.nd4j.linalg.api.ndarray.INDArray;

import java.util.List;

*/
/**
 * Encode a genomic position. Encode both chromosome and position with one hot encoding.
 * Created by fac2003 on 7/12/16.
 *//*

public class GenomicContextMapper implements FeatureMapper {
    int contextSize;
    private static final int NUM_POSITION_BITS = 32;
    private static final int NUM_CHROMOSOMES = 100;
    private final List<BinaryFeatureMapper> refContext;
    private ConcatFeatureMapper delegate;

    public GenomicContextMapper() {
        refContext = new ObjectArrayList<BinaryFeatureMapper>(21);
    }

    @Override
    public int numberOfFeatures() {
        return chromosomeMapper.numberOfFeatures() + positionMapper.numberOfFeatures();

    }

    @Override
    public void prepareToNormalize(BaseInformationRecords.BaseInformationOrBuilder record, int indexOfRecord) {
        contextSize = record.getGenomicSequenceContext().length();
        for (int i = 0; i < contextSize; i++){
            refContext.add(new BinaryFeatureMapper() {
                @Override
                public int getIntegerValue(BaseInformationRecords.BaseInformationOrBuilder record) {
                    return record.;
                }
            })
        }
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
*/
