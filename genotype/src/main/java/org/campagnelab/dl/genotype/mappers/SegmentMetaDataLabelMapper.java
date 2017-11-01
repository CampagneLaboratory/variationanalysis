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
 * Encode metadata about bases. One float per base, 0 when SNP matches reference, 1 when base has a snp, 2 when indel.
 * @author Fabien Campagne 11/1/17.
 */
public class SegmentMetaDataLabelMapper implements LabelMapper<SegmentInformationRecords.SegmentInformation>,
        ConfigurableLabelMapper {
    private int numberOfLabelsPerBase=1;
    private int maxSequenceLength;
    private int sampleIndex=0;

    @Override
    public int numberOfLabels() {
        return maxSequenceLength*numberOfLabelsPerBase;
    }

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
        return 0;
    }

    @Override
    public boolean hasMask() {
        return false;
    }

    @Override
    public void maskLabels(SegmentInformationRecords.SegmentInformation record, INDArray mask, int indexOfRecord) {

    }

    @Override
    public boolean isMasked(SegmentInformationRecords.SegmentInformation record, int featureIndex) {
        return true;
    }

    float[][] data;

    @Override
    public void prepareToNormalize(SegmentInformationRecords.SegmentInformation record, int indexOfRecord) {
        final int length = record.getLength();
        if (data == null) {
            data = new float[numberOfLabelsPerBase][maxSequenceLength];
        } else {
            for (int floatPerBaseIndex = 0; floatPerBaseIndex < numberOfLabelsPerBase; floatPerBaseIndex++) {
                Arrays.fill(data[floatPerBaseIndex], 0);
            }
        }
        // mark EOS as present by default:
        for (int baseIndex=length;baseIndex<maxSequenceLength;baseIndex++) {
            data[0][baseIndex] = 1;
        }
        SegmentInformationRecords.Sample sample = record.getSample(sampleIndex);
        int baseCount = Math.min(maxSequenceLength, sample.getBaseCount());

        for (int baseIndex = 0; baseIndex < baseCount; baseIndex++) {
            SegmentInformationRecords.Base base = sample.getBase(baseIndex);
            // NB: one less label is stored in protobuf than used in the model. The extra one in the model is for EOS

            for (int j = 0; j < numberOfLabelsPerBase; j++) {
                // NB: one metadata float per base, 0 when no variant, 1 when snp, 2 when indel.
                float metadataEncoding=0;
                metadataEncoding+=base.getIsVariant()?1:0;
                metadataEncoding+=base.getHasTrueIndel()?1:0;
                data[j][baseIndex] = metadataEncoding;
            }
        }
    }

    @Override
    public void configure(Properties readerProperties) {

        String maxSequenceLengthString = readerProperties.getProperty("genotypes.segments.maxSequenceLength");

        this.maxSequenceLength = Integer.parseInt(maxSequenceLengthString);
    }
}
