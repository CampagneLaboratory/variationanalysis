package org.campagnelab.dl.varanalysis.storage;

import org.apache.commons.io.FileUtils;
import org.campagnelab.dl.varanalysis.protobuf.BaseInformationRecords;
import org.junit.After;
import org.junit.Before;
import org.junit.Test;

import java.io.File;
import java.util.AbstractList;
import java.util.ArrayList;
import java.util.List;

import static org.junit.Assert.*;

/**
 * Created by mas2182 on 5/24/16.
 */
public class RecordWriterTest {

    private RecordWriter writer;
    private static String  filename = "test-results/genotypes_mutated_protofbuf.parquet";
    private static int blockSize = 256 * 1024 * 1024;
    private static int pageSize = 64 * 1024;


    @Before
    public void setUp() throws Exception {
        FileUtils.deleteQuietly(new File("test-results"));
        FileUtils.forceMkdir(new File("test-results"));
        writer = new RecordWriter(filename, blockSize, pageSize);
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
            builder.setReferenceIndex( i + 50);
            builder.setMutatedBase("mutatedBase");
            BaseInformationRecords.SampleCountInfo.Builder builderInfo = BaseInformationRecords.SampleCountInfo.newBuilder();
            builderInfo.setFromSequence("from");
            builderInfo.setToSequence("to");
            builderInfo.setMatchesReference(true);
            builderInfo.setGenotypeCountForwardStrand(1);
            builderInfo.setGenotypeCountReverseStrand(2);
            builder.addCounts(builderInfo.build());
            writer.writeRecord(builder.build());
        }

    }



}