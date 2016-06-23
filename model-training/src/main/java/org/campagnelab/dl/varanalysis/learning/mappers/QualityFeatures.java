package org.campagnelab.dl.varanalysis.learning.mappers;

import org.campagnelab.dl.model.utils.FeatureMapper;
import org.campagnelab.dl.varanalysis.learning.genotypes.BaseGenotypeCountFactory;
import org.campagnelab.dl.varanalysis.learning.genotypes.GenotypeCountFactory;
import org.campagnelab.dl.varanalysis.learning.iterators.AbstractFeatureMapper;
import org.campagnelab.dl.varanalysis.protobuf.BaseInformationRecords;
import org.campagnelab.dl.varanalysis.storage.RecordReader;
import org.nd4j.linalg.api.ndarray.INDArray;

import java.util.List;

/**
 * This is the quality score feature mapper. maps phred scores from parquet
 *
 * @author Remi Torracinta, rct66
 */

public class QualityFeatures extends AbstractFeatureMapper implements FeatureMapper {


    public static final int QUALITY_NORM = 1;


    public int numberOfFeatures() {
        // we need features for the normal sample and for the tumor sample:
        // multiply by additional 2 for sorted and unsorted features
        return MAX_GENOTYPES * 2 * 2 * 2;
    }

    public void prepareToNormalize(BaseInformationRecords.BaseInformationOrBuilder record, int indexOfRecord) {
        //shouldn't need to do anything
    }


    int[] indices = new int[]{0, 0};

    public void mapFeatures(BaseInformationRecords.BaseInformationOrBuilder record, INDArray inputs, int indexOfRecord) {
        indices[0] = indexOfRecord;
        for (int featureIndex = 0; featureIndex < numberOfFeatures(); featureIndex++) {
            indices[1] = featureIndex;
            inputs.putScalar(indices, produceFeature(record, featureIndex));
        }
    }

    public float produceFeature(BaseInformationRecords.BaseInformationOrBuilder record, int featureIndex) {
        return normalize(produceFeatureInternal(record, featureIndex), QUALITY_NORM);
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
        assert featureIndex >= 0 && featureIndex < MAX_GENOTYPES * 2 * 2 * 2 : "Only MAX_GENOTYPES*2*2 features";
        boolean sort = (featureIndex >= (MAX_GENOTYPES * 2 * 2));
        if (sort) featureIndex = featureIndex - (MAX_GENOTYPES * 2 * 2);
        if (featureIndex < MAX_GENOTYPES * 2) {
            // germline counts written first:
            final QualityGenotypeCount genotypeCount = (QualityGenotypeCount) getAllCounts(record, false, sort).get(featureIndex / 2);
            if ((featureIndex % 2) == 1) {
                // odd featureIndices are forward strand:
                return genotypeCount.getQualityScoreForward();
            } else {
                return genotypeCount.getQualityScoreReverse();
            }
        } else {
            // tumor counts written next:
            featureIndex -= MAX_GENOTYPES * 2;
            final QualityGenotypeCount genotypeCount = (QualityGenotypeCount) getAllCounts(record, true, sort).get(featureIndex / 2);
            if ((featureIndex % 2) == 1) {
                // odd featureIndices are forward strand:
                return genotypeCount.getQualityScoreForward();
            } else {
                return genotypeCount.getQualityScoreReverse();
            }
        }
    }

    @Override
    protected GenotypeCountFactory getGenotypeCountFactory() {
        return new BaseGenotypeCountFactory() {

            @Override
            public GenotypeCount create() {
                return new QualityGenotypeCount();
            }
        };
    }


    public static float avgQuality(List<Integer> list) {
        double sum = 0;
        for (Integer i : list)
            sum += Math.pow((double) 10, -((double) i / (double) 10));
        if (list.size() == 0) return 1;
        return (float) (sum / (double) list.size());
    }

    @Override
    protected void initializeCount(BaseInformationRecords.CountInfo sampleCounts, GenotypeCount count) {
        QualityGenotypeCount myCount = (QualityGenotypeCount) count;
        myCount.set(avgQuality(RecordReader.expandFreq(sampleCounts.getQualityScoresForwardStrandList())),
                avgQuality(RecordReader.expandFreq(sampleCounts.getQualityScoresReverseStrandList())));
    }
}

