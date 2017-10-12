package org.campagnelab.dl.genotype.mappers;

import org.campagnelab.dl.framework.mappers.*;
import org.campagnelab.dl.varanalysis.protobuf.SegmentInformationRecords;
import org.nd4j.linalg.api.ndarray.INDArray;

import java.util.Properties;

/**
 * Simply reads features from the .ssi and expose to DL4J as 2D tensor with mask.
 */
public class SingleBaseFeatureMapperV1 implements FeatureMapper<SegmentInformationRecords.SegmentInformation>, ConfigurableFeatureMapper {
    private final int sampleIndex;

    public SingleBaseFeatureMapperV1(int sampleIndex) {
        this.sampleIndex = sampleIndex;

    }

    @Override
    public int numberOfFeatures() {
        return numberOfFeaturesPerBase; // determine by the mapper used to produce the .ssi
    }

    private int maxSequenceLength = -1;
    private int numberOfFeaturesPerBase = -1;

    @Override
    public MappedDimensions dimensions() {
        return new MappedDimensions(maxSequenceLength, numberOfFeaturesPerBase);
    }

    int[] indices = new int[]{0, 0, 0};

    @Override
    public void mapFeatures(SegmentInformationRecords.SegmentInformation record, INDArray labels, int indexOfRecord) {
        assert maxSequenceLength > 0 : "maxSequenceLength must be positive";
        assert numberOfFeaturesPerBase > 0 : "numberOfFeaturesPerBase must be positive";

        indices[0] = indexOfRecord;
        for (int k = 0; k < numberOfFeatures(); k++) {
            indices[1] = k / numberOfFeaturesPerBase;
            indices[2] = k % numberOfFeaturesPerBase;
            labels.putScalar(indices, produceFeature(record, k));
        }

    }

    @Override
    public float produceFeature(SegmentInformationRecords.SegmentInformation record, int labelIndex) {
        return data[labelIndex / numberOfFeaturesPerBase][labelIndex % numberOfFeaturesPerBase];
    }


    @Override
    public boolean hasMask() {
        return true;
    }

    @Override
    public void maskFeatures(SegmentInformationRecords.SegmentInformation record, INDArray mask, int indexOfRecord) {
        for (int i = 0; i < record.getLength(); i++) {
            mask.putScalar(i, 1.0);
        }
    }

    @Override
    public boolean isMasked(SegmentInformationRecords.SegmentInformation record, int featureIndex) {
        int i = featureIndex / numberOfFeaturesPerBase;
        return i <= record.getLength();
    }

    float[][] data;

    @Override
    public void prepareToNormalize(SegmentInformationRecords.SegmentInformation record, int indexOfRecord) {

        final int length = record.getLength();
        if (data == null) {
            data = new float[length][numberOfFeaturesPerBase];
        }
        SegmentInformationRecords.Sample sample = record.getSample(sampleIndex);
        for (int i = 0; i < length; i++) {
            SegmentInformationRecords.Base base = sample.getBase(i);
            assert base.getFeaturesCount() == numberOfFeaturesPerBase :
                    String.format(
                            "the number of features per base must match between protobuf content and " +
                                    "genotypes.segments.numFeaturesPerBase property (in .ssip). " +
                                    "Found %d at base index=%d", base.getFeaturesCount(),
                            i);

            for (int j = 0; j < numberOfFeaturesPerBase; j++) {
                data[i][j] = base.getFeatures(j);
            }
        }
    }

    @Override
    public void configure(Properties readerProperties) {
        String nfpbString = readerProperties.getProperty("genotypes.segments.numFeaturesPerBase");

        numberOfFeaturesPerBase = Integer.parseInt(nfpbString);
        String maxSequenceLengthString = readerProperties.getProperty("genotypes.segments.maxSequenceLength");

        this.maxSequenceLength = Integer.parseInt(maxSequenceLengthString);
        System.out.printf("numFeaturesPerBase=%d, maxSequenceLength=%d %n", numberOfFeaturesPerBase, maxSequenceLength);


    }
}
