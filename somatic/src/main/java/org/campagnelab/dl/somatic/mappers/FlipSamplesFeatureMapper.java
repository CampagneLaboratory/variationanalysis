package org.campagnelab.dl.somatic.mappers;

import org.campagnelab.dl.framework.mappers.FeatureMapper;
import org.campagnelab.dl.framework.mappers.MappedDimensions;
import org.campagnelab.dl.varanalysis.protobuf.BaseInformationRecords;
import org.nd4j.linalg.api.ndarray.INDArray;

/**
 * Flip the samples of a record on the fly. Sample 1 becomes sample 0 and reciprocally. Tumor status is also switched.
 * Could be useful before calling mutator, but untested.
 * Created by fac2003 on 11/17/16.
 */
public class FlipSamplesFeatureMapper implements FeatureMapper<BaseInformationRecords.BaseInformation> {
    private FeatureMapper<BaseInformationRecords.BaseInformation> delegate;
    @Override
    public int numberOfFeatures() {
        return delegate.numberOfFeatures();
    }

    @Override
    public MappedDimensions dimensions() {
        return delegate.dimensions();
    }

    public void prepareToNormalize(BaseInformationRecords.BaseInformation record, int indexOfRecord) {

        delegate.prepareToNormalize(flip(record), indexOfRecord);
    }

    private BaseInformationRecords.BaseInformation flip(BaseInformationRecords.BaseInformation record) {
        assert record.getSamplesCount()==2:"Flip can only flip records with two samples.";
        BaseInformationRecords.BaseInformation.Builder builder = record.toBuilder();
        BaseInformationRecords.SampleInfo.Builder firstSample = builder.getSamples(0).toBuilder();
        BaseInformationRecords.SampleInfo.Builder secondSample = builder.getSamples(1).toBuilder();
        // flip isTumor status:
        firstSample.setIsTumor(!firstSample.getIsTumor());
        secondSample.setIsTumor(!secondSample.getIsTumor());
        builder.setSamples(0,secondSample);
        builder.setSamples(1,firstSample);
        return builder.build();
    }

    public void mapFeatures(BaseInformationRecords.BaseInformation record, INDArray inputs, int indexOfRecord) {
        delegate.mapFeatures(flip(record), inputs, indexOfRecord);
    }

    @Override
    public boolean hasMask() {
        return delegate.hasMask();
    }

    public void maskFeatures( BaseInformationRecords.BaseInformation record, INDArray mask, int indexOfRecord) {
       if (delegate.hasMask()) {
           delegate.maskFeatures(flip(record), mask, indexOfRecord);
       }
    }

    public boolean isMasked( BaseInformationRecords.BaseInformation record, int featureIndex) {
        if (delegate.hasMask()) {
            return delegate.isMasked(flip(record), featureIndex);
        }else {
            return false;
        }
    }

    public float produceFeature(BaseInformationRecords.BaseInformation record, int featureIndex) {
        return delegate.produceFeature(flip(record), featureIndex);
    }



}