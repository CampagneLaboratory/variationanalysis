package org.campagnelab.dl.varanalysis.storage;

import org.apache.commons.io.FileUtils;
import org.campagnelab.dl.varanalysis.protobuf.BaseInformationRecords;
import org.junit.After;
import org.junit.Before;
import org.junit.Test;

import java.io.File;


/**
 * Created by mas2182 on 5/24/16.
 */
public class RecordWriterTest {

    private RecordWriter writer;
    private static String filename = "test-results/genotypes_mutated_protobuf.parquet";
    private static int blockSize = 256 * 1024 * 1024;
    private static int pageSize = 64 * 1024;


    @Before
    public void setUp() throws Exception {
        FileUtils.deleteQuietly(new File("test-results"));
        FileUtils.forceMkdir(new File("test-results"));
        writer = new RecordWriter(filename, blockSize, pageSize, true);
    }

    @After
    public void tearDown() throws Exception {
        writer.close();
    }

    @Test
    public void writeRecords() throws Exception {

        for (int i = 1; i <= 100; i++) {
            BaseInformationRecords.BaseInformation.Builder builder = BaseInformationRecords.BaseInformation.newBuilder();
            builder.setMutated((i % 2 == 0));
            builder.setIndexOfMutatedBase(i + 100);
            builder.setPosition(i);
            builder.setReferenceBase("refbase");
            builder.setReferenceIndex(i + 50);
            builder.setMutatedBase("mutatedBase");

            //germline counts
            BaseInformationRecords.SampleInfo.Builder sampleBuilder = BaseInformationRecords.SampleInfo.newBuilder();
            BaseInformationRecords.CountInfo.Builder builderInfo = BaseInformationRecords.CountInfo.newBuilder();
            builderInfo.setFromSequence("from");
            builderInfo.setToSequence("to");
            builderInfo.setMatchesReference(true);
            builderInfo.setGenotypeCountForwardStrand(1);
            builderInfo.setGenotypeCountReverseStrand(2);
            sampleBuilder.addCounts(builderInfo.build());
            builder.addSamples(sampleBuilder.build());

            //somatic counts
            BaseInformationRecords.SampleInfo.Builder sampleBuilderS = BaseInformationRecords.SampleInfo.newBuilder();
            BaseInformationRecords.CountInfo.Builder builderInfoS = BaseInformationRecords.CountInfo.newBuilder();
            builderInfoS.setFromSequence("from");
            builderInfoS.setToSequence("to");
            builderInfoS.setMatchesReference(true);
            builderInfoS.setGenotypeCountForwardStrand(5);
            builderInfoS.setGenotypeCountReverseStrand(7);
            sampleBuilderS.addCounts(builderInfoS.build());
            builder.addSamples(sampleBuilderS.build());

            writer.writeRecord( builder.build());
        }

    }


}