package org.campagnelab.dl.framework.mappers.pretraining;

import org.campagnelab.dl.framework.mappers.ConcatFeatureMapper;
import org.campagnelab.dl.framework.mappers.TwoDimensionalConcatFeatureMapper;
import org.campagnelab.dl.framework.mappers.FeatureMapper;
import org.campagnelab.dl.framework.mappers.MappedDimensions;
import org.campagnelab.dl.framework.mappers.NAryFeatureMapper;
import org.campagnelab.dl.framework.mappers.RNNFeatureMapper;
import org.campagnelab.dl.framework.mappers.processing.TwoDimensionalRemoveMaskFeatureMapper;
import org.nd4j.linalg.api.ndarray.INDArray;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.util.function.Function;

/**
 * Concatenating + mask-removing feature mapper used for pretraining.
 *
 * Example mapping:
 * Sequence1:    ACTGTC
 * Sequence2:    ACTG
 * Sequence3:    AC
 *
 * Sequence1Map: ACTGTC!ACTGTC
 * Sequence2Map: ACTG!ACTG
 * Sequence3Map: AC!AC
 *
 * where "!" is an EOS token
 * I.e., maps sequences twice, with an EOS token in between
 *
 * Created by joshuacohen on 12/14/16.
 */
public class RNNPretrainingFeatureMapper<RecordType> implements FeatureMapper<RecordType> {
    static private Logger LOG = LoggerFactory.getLogger(RNNPretrainingFeatureMapper.class);

    public FeatureMapper<RecordType> delegate;

    public RNNPretrainingFeatureMapper(FeatureMapper<RecordType> domainMapper, Integer eosIndex,
                                       Function<RecordType, Integer> recordToSequenceLength) {
        MappedDimensions domainDimensions = domainMapper.dimensions();
        if (domainDimensions.numDimensions() != 2) {
            throw new IllegalArgumentException("Mapper must map two dimensional features");
        }
        int featuresPerTimeStep = domainDimensions.numElements(1);
        if (eosIndex != null && eosIndex > featuresPerTimeStep) {
            throw new IllegalArgumentException(String.format("Invalid EOS index %d greater than number of features %d",
                    eosIndex, featuresPerTimeStep));
        }
        int sequenceMapperPadding = (eosIndex != null && eosIndex == featuresPerTimeStep) || eosIndex == null ? 1 : 0;
        FeatureMapper<RecordType> sequenceMapper = new TwoDimensionalConcatFeatureMapper<>(sequenceMapperPadding,
                domainMapper);
        FeatureMapper<RecordType> eosMapper = createEosMapper(eosIndex,
                featuresPerTimeStep, recordToSequenceLength, LOG);
        delegate = new TwoDimensionalRemoveMaskFeatureMapper<>(new TwoDimensionalConcatFeatureMapper<>(sequenceMapper, sequenceMapper));
//        delegate = new TwoDimensionalConcatFeatureMapper<>(sequenceMapper, sequenceMapper);
//        delegate = new TwoDimensionalConcatFeatureMapper<>(sequenceMapper, eosMapper, sequenceMapper);
    }

    @Override
    public int numberOfFeatures() {
        return delegate.numberOfFeatures();
    }

    @Override
    public MappedDimensions dimensions() {
        return delegate.dimensions();
    }

    @Override
    public void prepareToNormalize(RecordType record, int indexOfRecord) {
        delegate.prepareToNormalize(record, indexOfRecord);
    }

    @Override
    public void mapFeatures(RecordType record, INDArray inputs, int indexOfRecord) {
        delegate.mapFeatures(record, inputs, indexOfRecord);
    }

    @Override
    public boolean hasMask() {
        return delegate.hasMask();
    }

    @Override
    public void maskFeatures(RecordType record, INDArray mask, int indexOfRecord) {
        delegate.maskFeatures(record, mask, indexOfRecord);
    }

    @Override
    public boolean isMasked(RecordType record, int featureIndex) {
        return delegate.isMasked(record, featureIndex);
    }

    @Override
    public float produceFeature(RecordType record, int featureIndex) {
        return delegate.produceFeature(record, featureIndex);
    }

    static <RecordType> FeatureMapper<RecordType> createEosMapper(Integer eosIndex, int featuresPerTimeStep,
                                                      Function<RecordType, Integer> recordToSequenceLength,
                                                      Logger log) {
        FeatureMapper<RecordType> eosConcatMapper;
        if (eosIndex != null) {
            int eosFeaturesStart = eosIndex;
            int eosFeaturesLeft = featuresPerTimeStep - eosIndex - 1;
            FeatureMapper<RecordType> eosMapperAtIndex = new NAryFeatureMapper<>(1,
                    false, true, true);
            FeatureMapper<RecordType> eosMapperBeforeIndex = new NAryFeatureMapper<>(eosFeaturesStart,
                    true, true, true);
            FeatureMapper<RecordType> eosMapperAfterIndex = new NAryFeatureMapper<>(eosFeaturesLeft,
                    true, true, true);
            if (eosFeaturesStart > 0 && eosFeaturesLeft > 0) {
                eosConcatMapper = new ConcatFeatureMapper<>(eosMapperBeforeIndex, eosMapperAtIndex, eosMapperAfterIndex);
            } else if (eosFeaturesStart > 0 && eosFeaturesLeft <= 0) {
                eosConcatMapper = new ConcatFeatureMapper<>(eosMapperBeforeIndex, eosMapperAtIndex);
            } else if (eosFeaturesStart <= 0 && eosFeaturesLeft > 0) {
                eosConcatMapper = new ConcatFeatureMapper<>(eosMapperAtIndex, eosMapperAfterIndex);
            } else {
                log.warn("EOS index may be invalid; no features before or after");
                eosConcatMapper = new ConcatFeatureMapper<>(eosMapperAtIndex);
            }
        } else {
            FeatureMapper<RecordType> eosMapperBeforeIndex = new NAryFeatureMapper<>(featuresPerTimeStep,
                    true, true, true);
            FeatureMapper<RecordType> eosMapperAtIndex = new NAryFeatureMapper<>(1,
                    false, true, true);
            eosConcatMapper = new ConcatFeatureMapper<>(eosMapperBeforeIndex, eosMapperAtIndex);
        }
        return new RNNFeatureMapper<>(recordToSequenceLength, eosConcatMapper);
    }
}
