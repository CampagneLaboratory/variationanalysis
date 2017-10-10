package org.campagnelab.dl.genotype.mappers;

import org.campagnelab.dl.framework.mappers.ConfigurableLabelMapper;
import org.campagnelab.dl.framework.mappers.LabelMapper;
import org.campagnelab.dl.framework.mappers.MappedDimensions;
import org.campagnelab.dl.varanalysis.protobuf.SegmentInformationRecords;
import org.nd4j.linalg.api.ndarray.INDArray;

import java.util.Properties;

/**
 *  Simply reads labels from the .ssi and expose to DL4J as 2D tensor with mask.
 */
public class SingleBaseLabelMapperV1 implements LabelMapper<SegmentInformationRecords.SegmentInformation>, ConfigurableLabelMapper {
    private final int sampleIndex;

    public SingleBaseLabelMapperV1(int sampleIndex) {
        this.sampleIndex = sampleIndex;
    }

    @Override
    public int numberOfLabels() {
        return numberOfLabelsPerBase; // Combinations of  A,C,T,G,-
    }

    private int ploidy = 2;
    private int maxSequenceLength = 300;
    private int numberOfLabelsPerBase = 15;

    @Override
    public MappedDimensions dimensions() {
        return new MappedDimensions(maxSequenceLength, numberOfLabels());
    }

    int[] indices = new int[]{0, 0, 0};

    @Override
    public void mapLabels(SegmentInformationRecords.SegmentInformation record, INDArray labels, int indexOfRecord) {
        indices[0] = indexOfRecord;
        for (int k = 0; k < numberOfLabels(); k++) {
            indices[1] = k / numberOfLabelsPerBase;
            indices[2] = k % numberOfLabelsPerBase;
            labels.putScalar(indices, produceLabel(record, k));
        }

    }

    @Override
    public float produceLabel(SegmentInformationRecords.SegmentInformation record, int labelIndex) {
        return data[labelIndex / numberOfLabelsPerBase][labelIndex % numberOfLabelsPerBase];
    }


    @Override
    public boolean hasMask() {
        return true;
    }

    @Override
    public void maskLabels(SegmentInformationRecords.SegmentInformation record, INDArray mask, int indexOfRecord) {
        for (int i = 0; i < record.getLength(); i++) {
            mask.putScalar(i, 1.0);
        }
    }

    @Override
    public boolean isMasked(SegmentInformationRecords.SegmentInformation record, int featureIndex) {
        int i = featureIndex / numberOfLabelsPerBase;
        int j = featureIndex % numberOfLabelsPerBase;
        return i <= record.getLength();
    }

    float[][] data;

    @Override
    public void prepareToNormalize(SegmentInformationRecords.SegmentInformation record, int indexOfRecord) {
        final int length = record.getLength();
        if (data == null) {
            data = new float[length][numberOfLabelsPerBase];
        }
        SegmentInformationRecords.Sample sample = record.getSample(sampleIndex);
        for (int i = 0; i < length; i++) {
            SegmentInformationRecords.Base base = sample.getBase(i);
            for (int j = 0; j < numberOfLabelsPerBase; j++) {
                data[i][j] = base.getLabels(j);
            }
        }
    }


    protected int getNumValues(int ploidy, int numAlleles) {
        if (ploidy == 0) {
            return numAlleles;
        }
        int increment = getNumValues(ploidy - 1, numAlleles);
        int count = 0;
        for (int i = 0; i < numAlleles; i++) {
            count += (--increment);
        }
        return count;
    }

    @Override
    public void configure(Properties readerProperties) {
        String ploidyString = readerProperties.getProperty("genotypes.ploidy");
        ploidy = Integer.parseInt(ploidyString);
        numberOfLabelsPerBase = getNumValues(ploidy, 5);
        String maxSequenceLengthString = readerProperties.getProperty("genotypes.segments.maxSequenceLength");

        this.maxSequenceLength = Integer.parseInt(maxSequenceLengthString);

        System.out.printf("Ploidy=%d, Number of alleles per base=%d maxSequenceLength=%d %n", ploidy,
                numberOfLabelsPerBase,maxSequenceLengthString);
    }
}
