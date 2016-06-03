package org.campagnelab.dl.varanalysis.learning.mappers;

import it.unimi.dsi.fastutil.objects.ObjectArrayList;
import org.campagnelab.dl.varanalysis.protobuf.BaseInformationRecords;
import org.nd4j.linalg.api.ndarray.INDArray;

import java.util.Collections;

/**
 * Estimates the number of distinct read indices per base position.
 * Created by fac2003 on 6/3/16.
 *
 * @author Fabien Campagne
 */
public class ReadIndexFeatures implements FeatureMapper {

    public static final int NUM_GENOTYPES = 5;
    public static final int NUM_SAMPLES = 2;

    @Override
    public int numberOfFeatures() {
        // only one feature per genotype: the number of distinct read indices in each sample. 5 genotypes max.
        return NUM_SAMPLES * NUM_GENOTYPES;
    }

    public static final int MAX_GENOTYPES = 5;

    int sumCounts;

    public void prepareToNormalize(BaseInformationRecords.BaseInformationOrBuilder record, int indexOfRecord) {
        indices[0] = indexOfRecord;
        sumCounts = 0;
        for (int featureIndex = 0; featureIndex < numberOfFeatures(); featureIndex++) {
            sumCounts += produceFeatureInternal(record, featureIndex);
        }
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
        return normalize(produceFeatureInternal(record, featureIndex), sumCounts);
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
            return getAllCounts(record, false).get(featureIndex / 2).getDistinctReadIndices();
        } else {
            // tumor counts written next:
            featureIndex -= MAX_GENOTYPES * 2;
            return getAllCounts(record, true).get(featureIndex / 2).getDistinctReadIndices();
        }
    }

    private boolean oneSampleHasTumor(java.util.List<org.campagnelab.dl.varanalysis.protobuf.BaseInformationRecords.SampleInfo> samples) {
        for (BaseInformationRecords.SampleInfo sample : samples) {
            if (sample.getIsTumor()) return true;
        }
        return false;

    }

    private ObjectArrayList<ReadIndexWithCounts> getAllCounts(BaseInformationRecords.BaseInformationOrBuilder record, boolean isTumor) {
        assert oneSampleHasTumor(record.getSamplesList()) : "at least one sample must have hasTumor=true.";
        ObjectArrayList<ReadIndexWithCounts> list = new ObjectArrayList();
        for (BaseInformationRecords.SampleInfo sampleInfo : record.getSamplesList()) {
            if (isTumor != sampleInfo.getIsTumor()) continue;
            for (BaseInformationRecords.CountInfo sampleCounts : sampleInfo.getCountsList()) {
                ReadIndexWithCounts count = new ReadIndexWithCounts(sampleCounts.getGenotypeCountForwardStrand(),
                        sampleCounts.getGenotypeCountReverseStrand(),
                        sampleCounts.getToSequence(),
                        sampleCounts.getReadIndicesForwardStrandList(),
                        sampleCounts.getReadIndicesReverseStrandList());
                list.add(count);
            }
        }
        // pad with zero until we have 10 elements:
        while (list.size() < MAX_GENOTYPES) {
            list.add(new ReadIndexWithCounts(0, 0, "N", Collections.EMPTY_LIST,Collections.EMPTY_LIST));
        }
        // trim the list at 5 elements because we will consider only the 5 genotypes with largest total counts:
        list.trim(MAX_GENOTYPES);
        //sort in decreasing order of counts:
        Collections.sort(list);
        return list;
    }

}
