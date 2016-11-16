package org.campagnelab.dl.somatic.learning.mappers.trio;

import org.campagnelab.dl.somatic.mappers.trio.MagnitudeFeatures2Trio;
import org.campagnelab.dl.varanalysis.protobuf.BaseInformationRecords;
import org.junit.Test;

import static org.junit.Assert.assertEquals;

/**
 * Created by fac2003 on 5/24/16.
 */
public class MagnitudeFeatures2TrioTest {
    @Test
    public void mapFeatures() throws Exception {
        MagnitudeFeatures2Trio magnitude = new MagnitudeFeatures2Trio();

        assertEquals(30, magnitude.numberOfFeatures());

        BaseInformationRecords.BaseInformation record = prepareRecord();
        magnitude.prepareToNormalize(record, 0);

        assertEquals(1f/6, magnitude.produceFeature(record, 0), 0.01);
        assertEquals(1f/11, magnitude.produceFeature(record, 1), 0.01);
        assertEquals(1f/1, magnitude.produceFeature(record, 2), 0.01);
        assertEquals(1f/4,magnitude.produceFeature(record, 10), 0.01);
        assertEquals(1f/18, magnitude.produceFeature(record, 11), 0.01);
        assertEquals(1f/1, magnitude.produceFeature(record, 12), 0.01);
        assertEquals(1f/5, magnitude.produceFeature(record, 20), 0.01);
        assertEquals(1f/3,magnitude.produceFeature(record, 21), 0.01);
        assertEquals(1f/1,magnitude.produceFeature(record, 22), 0.01);



    }


    private BaseInformationRecords.BaseInformation prepareRecord() {
        BaseInformationRecords.BaseInformation.Builder builder = BaseInformationRecords.BaseInformation.newBuilder();
        builder.setPosition(1);
        builder.setReferenceIndex(0);
        builder.setReferenceBase("A");


        //father counts
        BaseInformationRecords.SampleInfo.Builder sampleBuilder = BaseInformationRecords.SampleInfo.newBuilder();
        BaseInformationRecords.CountInfo.Builder builderInfo = BaseInformationRecords.CountInfo.newBuilder();
        builderInfo.setFromSequence("A");
        builderInfo.setToSequence("C");
        builderInfo.setMatchesReference(true);
        builderInfo.setGenotypeCountForwardStrand(10);
        builderInfo.setGenotypeCountReverseStrand(5);
        sampleBuilder.addCounts(builderInfo.build());
        builder.addSamples(sampleBuilder.build());

        //mother counts
        BaseInformationRecords.SampleInfo.Builder sampleBuilderM = BaseInformationRecords.SampleInfo.newBuilder();
        BaseInformationRecords.CountInfo.Builder builderInfoM = BaseInformationRecords.CountInfo.newBuilder();
        builderInfoM.setFromSequence("A");
        builderInfoM.setToSequence("C");
        builderInfoM.setMatchesReference(true);
        builderInfoM.setGenotypeCountForwardStrand(17);
        builderInfoM.setGenotypeCountReverseStrand(3);
        sampleBuilderM.addCounts(builderInfoM.build());
        builder.addSamples(sampleBuilderM.build());

        //somatic counts
        BaseInformationRecords.SampleInfo.Builder sampleBuilderS = BaseInformationRecords.SampleInfo.newBuilder();
        BaseInformationRecords.CountInfo.Builder builderInfoS = BaseInformationRecords.CountInfo.newBuilder();
        builderInfoS.setFromSequence("A");
        builderInfoS.setToSequence("C");
        builderInfoS.setMatchesReference(true);
        builderInfoS.setGenotypeCountForwardStrand(2);
        builderInfoS.setGenotypeCountReverseStrand(4);
        sampleBuilderS.addCounts(builderInfoS.build());
        sampleBuilderS.setIsTumor(true);
        builder.addSamples(sampleBuilderS.build());
        return builder.build();
    }

    //String[] expectedFeatures = {
     //       "[0.12, 0.12, 0.12, 0.12, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.12, 0.12, 0.06, 0.06, 0.00, 0.00, 0.06, 0.06, 0.00, 0.00]"};

}