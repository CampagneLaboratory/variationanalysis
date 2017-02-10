package org.campagnelab.dl.genotype.mappers;

import org.campagnelab.dl.framework.mappers.LabelMapper;
import org.campagnelab.dl.framework.mappers.MappedDimensions;
import org.campagnelab.dl.framework.mappers.RNNLabelMapper;
import org.campagnelab.dl.varanalysis.protobuf.BaseInformationRecords;
import org.nd4j.linalg.api.ndarray.INDArray;

/**
 * Created by joshuacohen on 2/8/17.
 */
public class TrueGenotypeLSTMLabelMapper implements LabelMapper<BaseInformationRecords.BaseInformation> {

    private final RNNLabelMapper<String> delegate;
    private static final int labelsPerTimeStep = 8;
    private String cachedRecordGenotype;

    public TrueGenotypeLSTMLabelMapper(int maxIndelLength) {
        delegate = new RNNLabelMapper<>(maxIndelLength, labelsPerTimeStep,
                TrueGenotypeLSTMLabelMapper::recordToLabel, String::length);
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
    public void mapLabels(BaseInformationRecords.BaseInformation record, INDArray labels, int indexOfRecord) {
        delegate.mapLabels(cachedRecordGenotype, labels, indexOfRecord);
    }

    @Override
    public float produceLabel(BaseInformationRecords.BaseInformation record, int labelIndex) {
        return delegate.produceLabel(cachedRecordGenotype, labelIndex);
    }

    @Override
    public boolean hasMask() {
        return delegate.hasMask();
    }

    @Override
    public void maskLabels(BaseInformationRecords.BaseInformation record, INDArray mask, int indexOfRecord) {
        delegate.maskLabels(cachedRecordGenotype, mask, indexOfRecord);
    }

    @Override
    public boolean isMasked(BaseInformationRecords.BaseInformation record, int featureIndex) {
        return delegate.isMasked(cachedRecordGenotype, featureIndex);
    }

    @Override
    public void prepareToNormalize(BaseInformationRecords.BaseInformation record, int indexOfRecord) {
        cachedRecordGenotype = record.getTrueGenotype();
    }

    private static int[] recordToLabel(String record) {
        return record.chars().map(TrueGenotypeLSTMLabelMapper::baseToLabel).toArray();
    }

    private static int baseToLabel(int base) {
        switch (base) {
            case 'a':
            case 'A':
                return 0;
            case 't':
            case 'T':
                return 1;
            case 'c':
            case 'C':
                return 2;
            case 'g':
            case 'G':
                return 3;
            case 'n':
            case 'N':
                return 4;
            case '-':
                return 5;
            case '/':
                return 6;
            default:
                return 7;

        }
    }
}
