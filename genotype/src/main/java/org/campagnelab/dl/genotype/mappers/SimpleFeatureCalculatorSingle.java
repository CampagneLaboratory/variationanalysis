package org.campagnelab.dl.genotype.mappers;

import org.campagnelab.dl.framework.mappers.FeatureCalculator;
import org.campagnelab.dl.somatic.genotypes.BaseGenotypeCountFactory;
import org.campagnelab.dl.somatic.genotypes.GenotypeCountFactory;
import org.campagnelab.dl.somatic.mappers.AbstractFeatureMapper;
import org.campagnelab.dl.somatic.mappers.GenotypeCount;
import org.campagnelab.dl.varanalysis.protobuf.BaseInformationRecords;
import org.nd4j.linalg.api.ndarray.INDArray;

/**
 * This is a simple feature mapper.
 * Each sample has the following information:
 * <pre> 69033640	11	false
 * position=14521   referenceIndex=0       isMutated=false
 * sample 0 counts=[10, 0, 0, 63, 0, 4, 0, 1, 62, 0]
 * <p>
 * Using these data, we can normalize counts and map them
 *
 * Created by rct66 on 11/16/16.
 *
 * @author Remi Torracinta
 */

public class SimpleFeatureCalculatorSingle extends AbstractFeatureMapper<BaseInformationRecords.BaseInformationOrBuilder>

        implements FeatureCalculator<BaseInformationRecords.BaseInformationOrBuilder> {


    public SimpleFeatureCalculatorSingle(boolean sort) {
        this.sort = sort;
    }

    public SimpleFeatureCalculatorSingle() {
        this.sort = true;
    }

    boolean normalized = false;


    @Override
    public int numberOfFeatures() {
        // we need features for the normal sample and for the tumor sample:

        return AbstractFeatureMapper.MAX_GENOTYPES * 2;
    }

    boolean sort; //by default we sort
    int sumCounts;

    public void prepareToNormalize(BaseInformationRecords.BaseInformationOrBuilder record, int indexOfRecord) {
        indices[0] = indexOfRecord;

        sumCounts = 0;
        for (int featureIndex = 0; featureIndex < numberOfFeatures(); featureIndex++) {
            sumCounts += produceFeatureInternal(record, featureIndex);
        }
        normalized = true;
    }

    @Override
    public int numberOfLabels() {
        return 2;
    }


    int[] indices = new int[]{0, 0};

    @Override
    public void mapFeatures(BaseInformationRecords.BaseInformationOrBuilder record, INDArray inputs, int indexOfRecord) {
        assert normalized : "prepareToNormalize must be called before mapFeatures.";
        indices[0] = indexOfRecord;

        final float[] buffer = getBuffer();
        mapFeatures(record, buffer, 0, indexOfRecord);
        for (int featureIndex = 0; featureIndex < numberOfFeatures(); featureIndex++) {
            indices[1] = featureIndex;
            inputs.putScalar(indices, buffer[featureIndex]);
        }
    }


    public void mapFeatures(BaseInformationRecords.BaseInformationOrBuilder record, float[] inputs, int offset, int indexOfRecord) {
        assert normalized : "prepareToNormalize must be called before mapFeatures.";
        for (int featureIndex = 0; featureIndex < numberOfFeatures(); featureIndex++) {
            inputs[featureIndex + offset] = produceFeature(record, featureIndex);
        }
    }

    public float produceFeature(BaseInformationRecords.BaseInformationOrBuilder record, int featureIndex) {
        return normalize(produceFeatureInternal(record, featureIndex), sumCounts);
    }

    @Override
    public String getFeatureName(int featureIndex) {
        assert (featureIndex >= 0 && featureIndex < AbstractFeatureMapper.MAX_GENOTYPES * 2) : "Only MAX_GENOTYPES*2*2 features";
        if ((featureIndex % 2) == 1) {
            // odd featureIndices are forward strand:
            return ("normGermlineForwardCount" + (featureIndex / 2));
        } else {
            return ("normGermlineReverseCount" + (featureIndex / 2));
        }
    }

    @Override
    public void mapLabels(BaseInformationRecords.BaseInformationOrBuilder record, INDArray labels, int indexOfRecord) {
        indices[0] = indexOfRecord;

        for (int labelIndex = 0; labelIndex < numberOfLabels(); labelIndex++) {
            indices[1] = labelIndex;
            labels.putScalar(indices, produceLabel(record, labelIndex));
        }
    }

    private float normalize(float value, int normalizationFactor) {
        if (normalizationFactor == 0) {
            return 0;
        }
        float normalized = value / normalizationFactor;
        assert normalized >= 0 && normalized <= 1 : "value must be normalized: " + normalized;
        return normalized;
    }


    public float produceFeatureInternal(BaseInformationRecords.BaseInformationOrBuilder record, int featureIndex) {
        assert (featureIndex >= 0 && featureIndex < AbstractFeatureMapper.MAX_GENOTYPES * 2) : "Only MAX_GENOTYPES*2*2 features";
        if ((featureIndex % 2) == 1) {
            // odd featureIndices are forward strand:
            return getAllCounts(record, false, sort).get(featureIndex / 2).forwardCount;
        } else {
            return getAllCounts(record, false, sort).get(featureIndex / 2).reverseCount;
        }
    }


    @Override
    public float produceLabel(BaseInformationRecords.BaseInformationOrBuilder record, int labelIndex) {
        assert labelIndex == 0 || labelIndex == 1 : "only one label.";

        //return record.getMutated() ? 1.0f : 0.0f;
        if (labelIndex == 0) return record.getMutated() ? 1 : 0;
        else {
            return !record.getMutated() ? 1 : 0;
        }

    }

    @Override
    protected void initializeCount(BaseInformationRecords.CountInfo sampleCounts, GenotypeCount count) {
        // nothing to do, already done in the base class.
    }

    @Override
    protected GenotypeCountFactory getGenotypeCountFactory() {

        return new BaseGenotypeCountFactory();

    }
}
