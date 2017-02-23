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
    static final int featuresOrLabelsPerTimeStep = 10;
    private String cachedRecordGenotype;
    private int maxGenotypeLength;

    public TrueGenotypeLSTMLabelMapper(int maxGenotypeLength) {
        this.maxGenotypeLength = maxGenotypeLength;
        delegate = new RNNLabelMapper<>(maxGenotypeLength + 2, featuresOrLabelsPerTimeStep,
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
        String trueGenotype = record.getTrueGenotype();
        StringBuilder cachedRecordGenotypeBuilder = new StringBuilder();
        cachedRecordGenotypeBuilder.append('$');
        if (trueGenotype.length() >= maxGenotypeLength) {
            cachedRecordGenotypeBuilder.append(trueGenotype.substring(0, maxGenotypeLength));
        } else {
            cachedRecordGenotypeBuilder.append(trueGenotype);
        }
        cachedRecordGenotypeBuilder.append('*');
        cachedRecordGenotype = cachedRecordGenotypeBuilder.toString();
        delegate.prepareToNormalize(cachedRecordGenotype, indexOfRecord);
    }

    private static int[] recordToLabel(String record) {
        return record.chars().map(TrueGenotypeLSTMLabelMapper::baseToLabel).toArray();
    }

    static int baseToLabel(int base) {
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
            case '|':
                return 6;
            case '*':
                return 7;
            case '$':
                return 8;
            default:
                return 9;

        }
    }
}
