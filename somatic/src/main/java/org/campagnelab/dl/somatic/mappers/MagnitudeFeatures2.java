package org.campagnelab.dl.somatic.mappers;

import org.campagnelab.dl.framework.mappers.FeatureCalculator;
import org.campagnelab.dl.somatic.genotypes.GenotypeCountFactory;
import org.campagnelab.dl.somatic.genotypes.BaseGenotypeCountFactory;
import org.campagnelab.dl.varanalysis.protobuf.BaseInformationRecords;
import org.nd4j.linalg.api.ndarray.INDArray;

/**
 *
 *
 * @author Fabien Campagne
 */

public class MagnitudeFeatures2 extends AbstractFeatureMapper<BaseInformationRecords.BaseInformationOrBuilder>
        implements FeatureCalculator<BaseInformationRecords.BaseInformationOrBuilder>
       {

    public MagnitudeFeatures2() {
    }

    @Override
    public int numberOfFeatures() {
        // we need features for the normal sample and for the tumor sample:

        return MAX_GENOTYPES * 2 * 2;
    }

    public void prepareToNormalize(BaseInformationRecords.BaseInformationOrBuilder record, int indexOfRecord) {
        indices[0] = indexOfRecord;
    }

    @Override
    public int numberOfLabels() {
        return 2;
    }

    int[] indices = new int[]{0, 0};

    @Override
    public void mapFeatures(BaseInformationRecords.BaseInformationOrBuilder record, INDArray inputs, int indexOfRecord) {
        indices[0] = indexOfRecord;
        prepareToNormalize(record, indexOfRecord);
        for (int featureIndex = 0; featureIndex < numberOfFeatures(); featureIndex++) {
            indices[1] = featureIndex;
            inputs.putScalar(indices, produceFeature(record, featureIndex));
        }
    }



           public float produceFeature(BaseInformationRecords.BaseInformationOrBuilder record, int featureIndex) {
        return normalize(produceFeatureInternal(record, featureIndex), 1);
    }

    @Override
    public String getFeatureName(int featureIndex) {
        assert (featureIndex >= 0 && featureIndex < MAX_GENOTYPES * 2 * 2) : "Only MAX_GENOTYPES*2*2 + 1 features";
        if (featureIndex < MAX_GENOTYPES * 2) {
            // germline counts written first:
            if ((featureIndex % 2) == 1) {
                // odd featureIndices are forward strand:
                return ("invGermlineForwardCount" + (featureIndex / 2));
            } else {
                return ("invGermlineReverseCount" + (featureIndex / 2));
            }
        } else {
            // tumor counts written next:
            featureIndex -= MAX_GENOTYPES * 2;
            if ((featureIndex % 2) == 1) {
                // odd featureIndices are forward strand:
                return ("invSomaticForwardCount" + (featureIndex / 2));
            } else {
                return ("invSomaticReverseCount" + (featureIndex / 2));
            }
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
        float normalized = 1 / (float) (value + 1);
        assert normalized >= 0 && normalized <= 1 : "value must be normalized: " + normalized;
        return normalized;
    }


    public float produceFeatureInternal(BaseInformationRecords.BaseInformationOrBuilder record, int featureIndex) {
        assert (featureIndex >= 0 && featureIndex < MAX_GENOTYPES * 2 * 2) : "Only MAX_GENOTYPES*2*2 + 1 features";
        if (featureIndex < MAX_GENOTYPES * 2) {
            // germline counts written first:
            if ((featureIndex % 2) == 1) {
                // odd featureIndices are forward strand:
                return getAllCounts(record, false, true).get(featureIndex / 2).forwardCount;
            } else {
                return getAllCounts(record, false, true).get(featureIndex / 2).reverseCount;
            }
        } else {
            // tumor counts written next:
            featureIndex -= MAX_GENOTYPES * 2;
            if ((featureIndex % 2) == 1) {
                // odd featureIndices are forward strand:
                return getAllCounts(record, true, true).get(featureIndex / 2).forwardCount;
            } else {
                return getAllCounts(record, true, true).get(featureIndex / 2).reverseCount;
            }
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
