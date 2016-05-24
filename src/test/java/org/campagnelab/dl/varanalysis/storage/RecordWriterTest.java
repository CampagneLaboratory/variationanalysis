package org.campagnelab.dl.varanalysis.storage;

import org.apache.commons.io.FileUtils;
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
    public void writeRecord() throws Exception {
        List<RecordWriter.SampleCountInfo> info = new ArrayList<RecordWriter.SampleCountInfo>();
        RecordWriter.SampleCountInfo info1 = new RecordWriter.SampleCountInfo();
        info1.fromSequence = "from";
        info1.toSequence = "to";
        info1.genotypeCountForwardStrand = 100;
        info1.genotypeCountReverseStrand = 200;
        info1.matchesReference = true;
        info.add(info1);
        RecordWriter.SampleCountInfo info2 = new RecordWriter.SampleCountInfo();
        info2.fromSequence = "from2";
        info2.toSequence = "to2";
        info2.genotypeCountForwardStrand = 300;
        info2.genotypeCountReverseStrand = 400;
        info2.matchesReference = false;
        info.add(info2);
        for (int i = 1; i <= 100; i++) {
            writer.writeRecord(info, "mutatedBase",
                    "refbase",  i + 50, i + 100, i, (i % 2 == 0));
        }

    }



}