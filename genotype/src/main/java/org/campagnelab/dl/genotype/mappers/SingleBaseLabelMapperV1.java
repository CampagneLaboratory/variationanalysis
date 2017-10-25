package org.campagnelab.dl.genotype.mappers;

import org.campagnelab.dl.framework.mappers.ConfigurableLabelMapper;
import org.campagnelab.dl.framework.mappers.LabelMapper;
import org.campagnelab.dl.framework.mappers.MappedDimensions;
import org.campagnelab.dl.varanalysis.protobuf.SegmentInformationRecords;
import org.nd4j.linalg.api.ndarray.INDArray;
import org.nd4j.linalg.factory.Nd4j;

import java.util.Arrays;
import java.util.Properties;

/**
 * Simply reads labels from the .ssi and expose to DL4J as 2D tensor with mask.
 */
public class SingleBaseLabelMapperV1 implements LabelMapper<SegmentInformationRecords.SegmentInformation>, ConfigurableLabelMapper {
    private final int sampleIndex;
    private float[] mask;

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


    @Override
    public void mapLabels(SegmentInformationRecords.SegmentInformation record, INDArray labels, int indexOfRecord) {
        INDArray dataCopy = Nd4j.create(data);
        dataCopy.detach();
        INDArray row = labels.tensorAlongDimension(indexOfRecord, 1, 2);
        row.assign(dataCopy);


    }

    @Override
    public float produceLabel(SegmentInformationRecords.SegmentInformation record, int labelIndex) {
        int colIndex = labelIndex % numberOfLabelsPerBase;
        int rowIndex = labelIndex / numberOfLabelsPerBase;
        return data[colIndex][rowIndex];
    }


    @Override
    public boolean hasMask() {
        return true;
    }

    @Override
    public void maskLabels(SegmentInformationRecords.SegmentInformation record, INDArray mask, int indexOfRecord) {
        INDArray recordMask=Nd4j.create(this.mask);
        recordMask.detach();
        INDArray row = mask.tensorAlongDimension(indexOfRecord, 1);
        row.assign(recordMask);

    }

    @Override
    public boolean isMasked(SegmentInformationRecords.SegmentInformation record, int featureIndex) {
        return this.mask[featureIndex]>0;
    }

    float[][] data;

    @Override
    public void prepareToNormalize(SegmentInformationRecords.SegmentInformation record, int indexOfRecord) {
        final int length = record.getLength();
        if (data == null) {
            data = new float[numberOfLabelsPerBase][maxSequenceLength];
            mask = new float[maxSequenceLength];
        } else {
            for (int floatPerBaseIndex = 0; floatPerBaseIndex < numberOfLabelsPerBase; floatPerBaseIndex++) {
                Arrays.fill(data[floatPerBaseIndex], 0);
            }
            Arrays.fill(mask,0);
        }
        SegmentInformationRecords.Sample sample = record.getSample(sampleIndex);
        int baseCount = Math.min(maxSequenceLength, sample.getBaseCount());

        for (int baseIndex = 0; baseIndex < baseCount; baseIndex++) {
            SegmentInformationRecords.Base base = sample.getBase(baseIndex);
            mask[baseIndex] = 1;
            assert base.getLabelsCount() == numberOfLabelsPerBase :
                    String.format(
                            "the number of labels per base must match between protobuf content and " +
                                    "genotypes.segments.numLabelsPerBase property (in .ssip). " +
                                    "Found %d at base index=%d", base.getLabelsCount(),
                            baseIndex);
            for (int j = 0; j < numberOfLabelsPerBase; j++) {
                data[j][baseIndex] = base.getLabels(j);
            }
        }
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
