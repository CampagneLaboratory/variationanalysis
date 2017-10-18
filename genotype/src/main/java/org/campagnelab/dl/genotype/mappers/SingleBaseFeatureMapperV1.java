package org.campagnelab.dl.genotype.mappers;

import org.campagnelab.dl.framework.mappers.*;
import org.campagnelab.dl.varanalysis.protobuf.SegmentInformationRecords;
import org.nd4j.linalg.api.ndarray.INDArray;
import scala.Array;

import java.util.Arrays;
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
        return numberOfFeaturesPerBase*maxSequenceLength; // determine by the mapper used to produce the .ssi
    }

    private int maxSequenceLength = -1;
    private int numberOfFeaturesPerBase = -1;

    @Override
    public MappedDimensions dimensions() {
        return new MappedDimensions( numberOfFeaturesPerBase,maxSequenceLength);
    }

    int[] indices = new int[]{0, 0, 0};

    @Override
    public void mapFeatures(SegmentInformationRecords.SegmentInformation record, INDArray labels, int indexOfRecord) {
        assert maxSequenceLength > 0 : "maxSequenceLength must be positive";
        assert numberOfFeaturesPerBase > 0 : "numberOfFeaturesPerBase must be positive";

        indices[0] = indexOfRecord;
        for (int k = 0; k < numberOfFeatures(); k++) {
            indices[1] = k % numberOfFeaturesPerBase;
            indices[2] = k / numberOfFeaturesPerBase;
            labels.putScalar(indices, produceFeature(record, k));
        }

    }

    @Override
    public float produceFeature(SegmentInformationRecords.SegmentInformation record, int featureIndex) {
        int colIndex = featureIndex % numberOfFeaturesPerBase;
        int rowIndex = featureIndex / numberOfFeaturesPerBase;
        return data[colIndex][rowIndex];
    }


    @Override
    public boolean hasMask() {
        return true;
    }

    @Override
    public void maskFeatures(SegmentInformationRecords.SegmentInformation record, INDArray mask, int indexOfRecord) {
        for (int i = 0; i < maxSequenceLength; i++) {
            mask.putScalar(indexOfRecord, i, (double)i<record.getLength()?1:0);
        }

    }

    @Override
    public boolean isMasked(SegmentInformationRecords.SegmentInformation record, int featureIndex) {
        int i = featureIndex / numberOfFeaturesPerBase;
        return i < record.getLength();
    }

    float[][] data;

    @Override
    public void prepareToNormalize(SegmentInformationRecords.SegmentInformation record, int indexOfRecord) {

        if (data == null) {
            data = new float[numberOfFeaturesPerBase][maxSequenceLength];
        }else {
            for (int baseIndex = 0; baseIndex < maxSequenceLength; baseIndex++) {
                for (int floatPerBaseIndex = 0; floatPerBaseIndex < numberOfFeaturesPerBase; floatPerBaseIndex++) {
                    data[floatPerBaseIndex][baseIndex] =0;
                }
            }
        }
        SegmentInformationRecords.Sample sample = record.getSample(sampleIndex);
        for (int baseIndex = 0; baseIndex < sample.getBaseCount(); baseIndex++) {
            SegmentInformationRecords.Base base = sample.getBase(baseIndex);
            assert base.getFeaturesCount() == numberOfFeaturesPerBase :
                    String.format(
                            "the number of features per base must match between protobuf content and " +
                                    "genotypes.segments.numFeaturesPerBase property (in .ssip). " +
                                    "Found %d at base index=%d", base.getFeaturesCount(),
                            baseIndex);

            for (int floatPerBaseIndex = 0; floatPerBaseIndex < numberOfFeaturesPerBase; floatPerBaseIndex++) {
                data[floatPerBaseIndex][baseIndex] = base.getFeatures(floatPerBaseIndex);
            }
        }
    }

    @Override
    public void configure(Properties readerProperties) {
        String nfpbString = readerProperties.getProperty("genotypes.segments.numFeaturesPerBase");
        assert nfpbString != null : "The .ssip file must define property: genotypes.segments.numFeaturesPerBase";

        numberOfFeaturesPerBase = Integer.parseInt(nfpbString);
        String maxSequenceLengthString = readerProperties.getProperty("genotypes.segments.maxSequenceLength");
        assert maxSequenceLengthString != null : "The .ssip file must define property: genotypes.segments.maxSequenceLength";

        this.maxSequenceLength = Integer.parseInt(maxSequenceLengthString);
        System.out.printf("numFeaturesPerBase=%d, maxSequenceLength=%d %n", numberOfFeaturesPerBase, maxSequenceLength);


    }
}
