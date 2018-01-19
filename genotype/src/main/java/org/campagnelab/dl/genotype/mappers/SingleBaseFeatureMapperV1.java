package org.campagnelab.dl.genotype.mappers;

import org.campagnelab.dl.framework.mappers.ConfigurableFeatureMapper;
import org.campagnelab.dl.framework.mappers.FeatureMapper;
import org.campagnelab.dl.framework.mappers.MappedDimensions;
import org.campagnelab.dl.varanalysis.protobuf.SegmentInformationRecords;
import org.nd4j.linalg.api.ndarray.INDArray;
import org.nd4j.linalg.factory.Nd4j;

import java.util.Arrays;
import java.util.Properties;

/**
 * Simply reads features from the .ssi and expose to DL4J as 2D tensor with mask.
 */
public class SingleBaseFeatureMapperV1 implements FeatureMapper<SegmentInformationRecords.SegmentInformation>, ConfigurableFeatureMapper
{
    private int sampleIndex;
    private float[] mask;

    public SingleBaseFeatureMapperV1(int sampleIndex) {
        this.sampleIndex = sampleIndex;
    }

    public SingleBaseFeatureMapperV1() {
        this(0);
    }

    @Override
    public void setSampleIndex(int sampleIndex) {
        this.sampleIndex=sampleIndex;
    }

    @Override
    public int numberOfFeatures() {
        return numberOfFeaturesPerBase * maxSequenceLength; // determine by the mapper used to produce the .ssi
    }

    private int maxSequenceLength = -1;
    private int numberOfFeaturesPerBase = -1;

    @Override
    public MappedDimensions dimensions() {
        return new MappedDimensions(numberOfFeaturesPerBase, maxSequenceLength);
    }

    int[] indices = new int[]{0, 0, 0};

    @Override
    public void mapFeatures(SegmentInformationRecords.SegmentInformation record, INDArray features, int indexOfRecord) {
        assert maxSequenceLength > 0 : "maxSequenceLength must be positive";
        assert numberOfFeaturesPerBase > 0 : "numberOfFeaturesPerBase must be positive";

        INDArray dataCopy = Nd4j.create(data);
        dataCopy.detach();
        INDArray row = features.tensorAlongDimension(indexOfRecord, 1, 2);
        row.assign(dataCopy);
  //      System.out.printf("indexOfRecord: %d features:\n%s",indexOfRecord,features);
    //    System.out.println("---");

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

        if (data == null) {
            data = new float[numberOfFeaturesPerBase][maxSequenceLength];
            mask = new float[maxSequenceLength];
        } else {
            for (int floatPerBaseIndex = 0; floatPerBaseIndex < numberOfFeaturesPerBase; floatPerBaseIndex++) {
                Arrays.fill(data[floatPerBaseIndex], 0);
            }
            Arrays.fill(mask,0);
        }
        SegmentInformationRecords.Sample sample = record.getSample(sampleIndex);
        final int numberOfBases = Math.min(maxSequenceLength, sample.getBaseCount());
        for (int baseIndex = 0; baseIndex < numberOfBases; baseIndex++) {
            SegmentInformationRecords.Base base = sample.getBase(baseIndex);
            mask[baseIndex] = 1;
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
