package org.campagnelab.dl.varanalysis.learning.iterators;

import org.campagnelab.dl.varanalysis.learning.mappers.MagnitudeFeatures;
import org.campagnelab.dl.varanalysis.protobuf.BaseInformationRecords;
import org.junit.Test;

import static org.junit.Assert.*;

/**
 * Created by fac2003 on 5/24/16.
 */
public class MagnitudeFeaturesTest {
    @Test
    public void mapFeatures() throws Exception {
        MagnitudeFeatures magnitude = new MagnitudeFeatures();

        assertEquals(31, magnitude.numberOfFeatures());

        BaseInformationRecords.BaseInformation record = prepareRecord(33);
        magnitude.prepareToNormalize(record, 0);

        assertEquals(1, magnitude.produceFeature(record, 0), 0.01);
        assertEquals(1, magnitude.produceFeature(record, 1), 0.01);
        assertEquals(0, magnitude.produceFeature(record, 2), 0.01);
        assertEquals(0, magnitude.produceFeature(record, 3), 0.01);

        record = prepareRecord(151);
        magnitude.prepareToNormalize(record, 0);

        assertEquals(1, magnitude.produceFeature(record, 0), 0.01);
        assertEquals(1, magnitude.produceFeature(record, 1), 0.01);
        assertEquals(1, magnitude.produceFeature(record, 2), 0.01);
        assertEquals(0, magnitude.produceFeature(record, 3), 0.01);


    }


    private BaseInformationRecords.BaseInformation prepareRecord(int sumTotal) {
        BaseInformationRecords.BaseInformation.Builder builder = BaseInformationRecords.BaseInformation.newBuilder();
        builder.setPosition(1);
        builder.setReferenceIndex(0);
        builder.setReferenceBase("A");


        //germline counts
        BaseInformationRecords.SampleInfo.Builder sampleBuilder = BaseInformationRecords.SampleInfo.newBuilder();
        BaseInformationRecords.CountInfo.Builder builderInfo = BaseInformationRecords.CountInfo.newBuilder();
        builderInfo.setFromSequence("A");
        builderInfo.setToSequence("C");
        builderInfo.setMatchesReference(true);
        builderInfo.setGenotypeCountForwardStrand(sumTotal/4);
        builderInfo.setGenotypeCountReverseStrand(sumTotal/4);
        sampleBuilder.addCounts(builderInfo.build());
        builder.addSamples(sampleBuilder.build());

        //somatic counts
        BaseInformationRecords.SampleInfo.Builder sampleBuilderS = BaseInformationRecords.SampleInfo.newBuilder();
        BaseInformationRecords.CountInfo.Builder builderInfoS = BaseInformationRecords.CountInfo.newBuilder();
        builderInfoS.setFromSequence("A");
        builderInfoS.setToSequence("T");
        builderInfoS.setMatchesReference(true);
        builderInfoS.setGenotypeCountForwardStrand(sumTotal/4);
        builderInfoS.setGenotypeCountReverseStrand(sumTotal/4);
        sampleBuilderS.addCounts(builderInfoS.build());
        builder.addSamples(sampleBuilderS.build());
        return builder.build();
    }

}