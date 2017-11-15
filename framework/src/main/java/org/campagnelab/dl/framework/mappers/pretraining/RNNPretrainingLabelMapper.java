package org.campagnelab.dl.framework.mappers.pretraining;

import org.campagnelab.dl.framework.mappers.*;
import org.campagnelab.dl.framework.mappers.processing.TwoDimensionalRemoveMaskFeatureMapper;
import org.nd4j.linalg.api.ndarray.INDArray;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.util.function.Function;

/**
 * Concatenating + mask-removing label mapper used for pretraining.
 *
 * Example mapping:
 * Sequence1:    ACTGTC
 * Sequence2:    ACTG
 * Sequence3:    AC
 *
 * Sequence1Map: ......ACTGTC!
 * Sequence2Map: ....ACTG!
 * Sequence3Map: ..AC!
 *
 * where "!" is an EOS token and "." is a blank character token
 * I.e., maps blank tokens of the length of the sequence, followed by the sequence, followed by an EOS token
 *
 * Created by joshuacohen on 12/14/16.
 */
public class RNNPretrainingLabelMapper<RecordType> implements LabelMapper<RecordType> {
    static private Logger LOG = LoggerFactory.getLogger(RNNPretrainingLabelMapper.class);

    public LabelMapper<RecordType> delegate;

    public RNNPretrainingLabelMapper(FeatureMapper<RecordType> domainMapper, Integer eosIndex,
                                       Function<RecordType, Integer> recordToSequenceLength) {
        MappedDimensions domainDimensions = domainMapper.dimensions();
        if (domainDimensions.numDimensions() != 2) {
            throw new IllegalArgumentException("Mapper must map two dimensional labels");
        }
        int numTimeSteps = domainDimensions.numElements(2);
        int labelsPerTimeStep = domainDimensions.numElements(1);
        if (eosIndex != null && eosIndex > labelsPerTimeStep) {
            throw new IllegalArgumentException(String.format("Invalid EOS index %d greater than number of features %d",
                    eosIndex, labelsPerTimeStep));
        }
        int sequenceMapperPadding = (eosIndex != null && eosIndex == labelsPerTimeStep) || eosIndex == null ? 1 : 0;
        FeatureMapper<RecordType>[] zeroMapperDelegates = new FeatureMapper[numTimeSteps];
        for (int i = 0; i < numTimeSteps; i++) {
            zeroMapperDelegates[i] = new NAryFeatureMapper<>(labelsPerTimeStep + sequenceMapperPadding,
                    true, true, false);
        }
        FeatureMapper<RecordType> sequenceMapper = new TwoDimensionalConcatFeatureMapper<>(sequenceMapperPadding,
                domainMapper);
        FeatureMapper<RecordType> zeroMapper = new RNNFeatureMapper<>(recordToSequenceLength, zeroMapperDelegates);
        FeatureMapper<RecordType> eosMapper = RNNPretrainingFeatureMapper.createEosMapper(eosIndex,
                labelsPerTimeStep, recordToSequenceLength, LOG);
        delegate = new LabelFromFeatureMapper<>(new TwoDimensionalRemoveMaskFeatureMapper<>(
                new TwoDimensionalConcatFeatureMapper<>(zeroMapper, sequenceMapper, eosMapper),
                recordToSequenceLength));
    }

    @Override
    public int numberOfLabels() {
        return delegate.numberOfLabels();
    }

    @Override
    public MappedDimensions dimensions() {
        return delegate.dimensions();
    }

    @Override
    public void mapLabels(RecordType record, INDArray labels, int indexOfRecord) {
        delegate.mapLabels(record, labels, indexOfRecord);
    }

    @Override
    public float produceLabel(RecordType record, int labelIndex) {
        return delegate.produceLabel(record, labelIndex);
    }

    @Override
    public boolean hasMask() {
        return delegate.hasMask();
    }

    @Override
    public void maskLabels(RecordType record, INDArray mask, int indexOfRecord) {
        delegate.maskLabels(record, mask, indexOfRecord);
    }

    @Override
    public boolean isMasked(RecordType record, int labelIndex) {
        return delegate.isMasked(record, labelIndex);
    }

    @Override
    public void prepareToNormalize(RecordType record, int indexOfRecord) {
        delegate.prepareToNormalize(record, indexOfRecord);
    }
}
