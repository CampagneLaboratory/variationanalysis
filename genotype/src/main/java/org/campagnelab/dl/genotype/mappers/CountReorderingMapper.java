package org.campagnelab.dl.genotype.mappers;

import org.campagnelab.dl.framework.mappers.AbstractFeatureMapper1D;
import org.campagnelab.dl.framework.mappers.FeatureNameMapper;
import org.campagnelab.dl.varanalysis.protobuf.BaseInformationRecords;

/**
 * Created by fac2003 on 12/15/16.
 */
public class CountReorderingMapper extends AbstractFeatureMapper1D<BaseInformationRecords.BaseInformationOrBuilder> {
    private int sampleIndex;
    private FeatureNameMapper<BaseInformationRecords.BaseInformationOrBuilder> delegate;

    public CountReorderingMapper(int sampleIndex, FeatureNameMapper<BaseInformationRecords.BaseInformationOrBuilder> delegate) {
        this.delegate = delegate;
        this.sampleIndex=sampleIndex;
    }

    @Override
    public String getFeatureName(int featureIndex) {
        return delegate.getFeatureName(featureIndex);
    }

    @Override
    public int numberOfFeatures() {
        return delegate.numberOfFeatures();
    }

    private BaseInformationRecords.BaseInformationOrBuilder sortedCountRecord;
    private RecordCountSortHelper sortHelper = new RecordCountSortHelper();

    @Override
    public void prepareToNormalize(BaseInformationRecords.BaseInformationOrBuilder record, int indexOfRecord) {

        sortedCountRecord = sortHelper.sort(sampleIndex,record);
        delegate.prepareToNormalize(sortedCountRecord, indexOfRecord);
    }

    @Override
    public float produceFeature(BaseInformationRecords.BaseInformationOrBuilder record, int featureIndex) {
        return delegate.produceFeature(sortedCountRecord, featureIndex);
    }
}
