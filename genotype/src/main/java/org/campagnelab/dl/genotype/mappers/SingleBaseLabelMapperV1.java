package org.campagnelab.dl.genotype.mappers;

import org.campagnelab.dl.framework.mappers.ConfigurableLabelMapper;
import org.campagnelab.dl.framework.mappers.LabelMapper;
import org.campagnelab.dl.framework.mappers.MappedDimensions;
import org.campagnelab.dl.varanalysis.protobuf.SegmentInformationRecords;
import org.nd4j.linalg.api.ndarray.INDArray;

import java.util.Properties;

/**
 * Simply reads labels from the .ssi and expose to DL4J as 2D tensor with mask.
 */
public class SingleBaseLabelMapperV1 implements LabelMapper<SegmentInformationRecords.SegmentInformation>, ConfigurableLabelMapper {
    private final int sampleIndex;

    public SingleBaseLabelMapperV1(int sampleIndex) {
        this.sampleIndex = sampleIndex;
    }

    @Override
    public int numberOfLabels() {
        return numberOfLabelsPerBase * maxSequenceLength; // Combinations of  A,C,T,G,-
    }

    private int ploidy = 2;
    private int maxSequenceLength = 300;
    private int numberOfLabelsPerBase = 15;

    @Override
    public MappedDimensions dimensions() {
        return new MappedDimensions(numberOfLabelsPerBase, maxSequenceLength);
    }

    int[] indices = new int[]{0, 0, 0};

    @Override
    public void mapLabels(SegmentInformationRecords.SegmentInformation record, INDArray labels, int indexOfRecord) {
        indices[0] = indexOfRecord;
        for (int k = 0; k < numberOfLabels(); k++) {
            indices[2] = k % numberOfLabelsPerBase;
            indices[1] = k / numberOfLabelsPerBase;
            labels.putScalar(indices, produceLabel(record, k));
        }

    }

    @Override
    public float produceLabel(SegmentInformationRecords.SegmentInformation record, int labelIndex) {
        return data[labelIndex % numberOfLabelsPerBase][labelIndex / numberOfLabelsPerBase];
    }


    @Override
    public boolean hasMask() {
        return true;
    }

    @Override
    public void maskLabels(SegmentInformationRecords.SegmentInformation record, INDArray mask, int indexOfRecord) {
        for (int i = 0; i < maxSequenceLength; i++) {
            mask.putScalar(indexOfRecord, i, (double) i < record.getLength() ? 1 : 0);
        }

    }

    @Override
    public boolean isMasked(SegmentInformationRecords.SegmentInformation record, int featureIndex) {
        int i = featureIndex / numberOfLabelsPerBase;
        int j = featureIndex % numberOfLabelsPerBase;
        return i < record.getLength();
    }

    float[][] data;

    @Override
    public void prepareToNormalize(SegmentInformationRecords.SegmentInformation record, int indexOfRecord) {
        final int length = record.getLength();
        if (data == null) {
            data = new float[numberOfLabelsPerBase][maxSequenceLength];
        } else {
            for (int baseIndex = 0; baseIndex < maxSequenceLength; baseIndex++) {
                for (int floatPerBaseIndex = 0; floatPerBaseIndex < numberOfLabelsPerBase; floatPerBaseIndex++) {
                    data[floatPerBaseIndex][baseIndex] = 0;
                }
            }
        }
        SegmentInformationRecords.Sample sample = record.getSample(sampleIndex);
        int baseCount = Math.min(maxSequenceLength, sample.getBaseCount());
        for (int i = 0; i < baseCount; i++) {
            SegmentInformationRecords.Base base = sample.getBase(i);
            for (int j = 0; j < numberOfLabelsPerBase; j++) {
                data[j][i] = base.getLabels(j);
            }
        }
    }


    public static int getNumValues(int ploidy, int numAlleles) {
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
        String maxNumberLabelsString = readerProperties.getProperty("maxNumOfLabels");

        this.numberOfLabelsPerBase = Integer.parseInt(maxNumberLabelsString);

        String maxSequenceLengthString = readerProperties.getProperty("genotypes.segments.maxSequenceLength");

        this.maxSequenceLength = Integer.parseInt(maxSequenceLengthString);

        System.out.printf("Ploidy=%d, Number of alleles per base=%d maxSequenceLength=%d %n", ploidy,
                numberOfLabelsPerBase, maxSequenceLength);
    }
}
