package org.campagnelab.dl.varanalysis.learning.mappers;

import it.unimi.dsi.fastutil.objects.ObjectArrayList;
import org.campagnelab.dl.varanalysis.learning.genotypes.BaseGenotypeCountFactory;
import org.campagnelab.dl.varanalysis.learning.genotypes.GenotypeCountFactory;
import org.campagnelab.dl.varanalysis.learning.iterators.AbstractFeatureMapper;
import org.campagnelab.dl.varanalysis.protobuf.BaseInformationRecords;
import org.campagnelab.dl.varanalysis.storage.RecordReader;
import org.nd4j.linalg.api.ndarray.INDArray;

import java.util.Collections;

/**
 * Estimates the number of distinct read indices per base position.
 * Created by fac2003 on 6/3/16.
 *
 * @author Fabien Campagne
 */
public class ReadIndexFeatures extends AbstractFeatureMapper implements FeatureMapper {


    @Override
    public int numberOfFeatures() {
        // only one feature per genotype: the number of distinct read indices in each sample. 5 genotypes max.
        // sorted and unsorted
        return MAX_GENOTYPES * 2 * 2 * 2;
    }

    int sumCounts;

    public void prepareToNormalize(BaseInformationRecords.BaseInformationOrBuilder record, int indexOfRecord) {
        indices[0] = indexOfRecord;
        sumCounts = 0;
        for (int featureIndex = 0; featureIndex < numberOfFeatures()/2; featureIndex++) {
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
        assert featureIndex >= 0 && featureIndex < MAX_GENOTYPES * 2 * 2 * 2 : "Only MAX_GENOTYPES*2*2 features";
        boolean sort = (featureIndex >= (MAX_GENOTYPES * 2 * 2));
        if (sort) featureIndex = featureIndex - (MAX_GENOTYPES * 2 * 2);
        if (featureIndex < MAX_GENOTYPES * 2) {
            // germline counts written first:
            final ReadIndexWithCounts genotypeCount = (ReadIndexWithCounts) getAllCounts(record, false,sort).get(featureIndex / 2);
            return genotypeCount.getDistinctReadIndices();
        } else {
            // tumor counts written next:
            featureIndex -= MAX_GENOTYPES * 2;
            final ReadIndexWithCounts genotypeCount = (ReadIndexWithCounts) getAllCounts(record, true,sort).get(featureIndex / 2);
            return genotypeCount.getDistinctReadIndices();
        }
    }

    private boolean oneSampleHasTumor(java.util.List<org.campagnelab.dl.varanalysis.protobuf.BaseInformationRecords.SampleInfo> samples) {
        for (BaseInformationRecords.SampleInfo sample : samples) {
            if (sample.getIsTumor()) return true;
        }
        return false;

    }

    @Override
    protected void initializeCount(BaseInformationRecords.CountInfo sampleCounts, GenotypeCount count) {
        ReadIndexWithCounts myCounts = (ReadIndexWithCounts) count;
        myCounts.set(RecordReader.expandFreq(sampleCounts.getReadIndicesForwardStrandList()),
                RecordReader.expandFreq(sampleCounts.getReadIndicesReverseStrandList()));
    }

    @Override
    protected GenotypeCountFactory getGenotypeCountFactory() {

        return new BaseGenotypeCountFactory() {
            @Override
            public GenotypeCount create() {
                return new ReadIndexWithCounts();
            }
        };
    }

}
