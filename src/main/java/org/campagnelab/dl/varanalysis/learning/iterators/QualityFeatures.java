package org.campagnelab.dl.varanalysis.learning.iterators;

import it.unimi.dsi.fastutil.objects.ObjectArrayList;
import org.campagnelab.dl.varanalysis.learning.mappers.FeatureMapper;
import org.campagnelab.dl.varanalysis.protobuf.BaseInformationRecords;
import org.nd4j.linalg.api.ndarray.INDArray;

import java.util.Collections;
import java.util.List;

/**
 * This is the quality score feature mapper. maps phred scores from parquet
 *
 * @author Remi Torracinta, rct66
 */

public class QualityFeatures implements FeatureMapper {


    public static final int MAX_GENOTYPES = 5;
    public static final int QUALITY_NORM = 1;


    public int numberOfFeatures() {
        // we need features for the normal sample and for the tumor sample:

        return MAX_GENOTYPES * 2*2;
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
        assert featureIndex >= 0 && featureIndex < MAX_GENOTYPES * 2*2 : "Only MAX_GENOTYPES*2*2 features";
        if (featureIndex < MAX_GENOTYPES * 2) {
            // germline counts written first:
            if ((featureIndex % 2) == 1) {
                // odd featureIndices are forward strand:
                return getAllCounts(record, false).get(featureIndex / 2).getQualityScoreForward();
            } else {
                return getAllCounts(record, false).get(featureIndex / 2).getQualityScoreReverse();
            }
        } else {
            // tumor counts written next:
            featureIndex -= MAX_GENOTYPES * 2;
            if ((featureIndex % 2) == 1) {
                // odd featureIndices are forward strand:
                return getAllCounts(record, true).get(featureIndex / 2).getQualityScoreForward();
            } else {
                return getAllCounts(record, true).get(featureIndex / 2).getQualityScoreReverse();
            }
        }
    }

    private ObjectArrayList<QualityGenotypeCount> getAllCounts(BaseInformationRecords.BaseInformationOrBuilder record, boolean isTumor) {
        assert oneSampleHasTumor(record.getSamplesList()) : "at least one sample must have hasTumor=true.";
        ObjectArrayList<QualityGenotypeCount> list = new ObjectArrayList();
        for (BaseInformationRecords.SampleInfo sampleInfo : record.getSamplesList()) {
            if (isTumor != sampleInfo.getIsTumor()) continue;
            for (BaseInformationRecords.CountInfo sampleCounts : sampleInfo.getCountsList()) {
                assert ((sampleCounts.getQualityScoresForwardStrandCount() > 0 ) || (sampleCounts.getQualityScoresReverseStrandCount() > 0))
                        : "record has no quality scores.";
                QualityGenotypeCount count = new QualityGenotypeCount(
                        sampleCounts.getGenotypeCountForwardStrand(),
                        sampleCounts.getGenotypeCountReverseStrand(),
                        sampleCounts.getToSequence(),
                        avgQuality(sampleCounts.getQualityScoresForwardStrandList()),
                        avgQuality(sampleCounts.getQualityScoresReverseStrandList()));
                list.add(count);
            }
        }
        // pad with zero until we have 10 elements:
        while (list.size() < MAX_GENOTYPES) {
            list.add(new QualityGenotypeCount(0, 0, "N", 0, 0));
        }
        // trim the list at 5 elements because we will consider only the 5 genotypes with largest total counts:
        list.trim(MAX_GENOTYPES);
        //sort in decreasing order of counts:
        Collections.sort(list);
        return list;
    }


    public static float avgQuality(List<Integer> list){
        double sum = 0;
        for (Integer i : list)
            sum += Math.pow((double)10, -((double)i/(double)10));
        return (float) (sum/(double)list.size());
    }


    private boolean oneSampleHasTumor(java.util.List<org.campagnelab.dl.varanalysis.protobuf.BaseInformationRecords.SampleInfo> samples) {
        for (BaseInformationRecords.SampleInfo sample : samples) {
            if (sample.getIsTumor()) return true;
        }
        return false;

    }
 }

