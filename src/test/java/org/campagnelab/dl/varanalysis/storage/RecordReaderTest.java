package org.campagnelab.dl.varanalysis.storage;

import org.junit.After;
import org.junit.Before;
import org.junit.Test;

import static org.junit.Assert.*;

/**
 * Created by mas2182 on 5/24/16.
 */
public class RecordReaderTest {
    private String  filename = "sample_data/genotypes_mutated.parquet";

    private RecordReader reader;

    @Before
    public void setUp() throws Exception {
         reader = new RecordReader(filename);
    }

    @After
    public void tearDown() throws Exception {
        reader.close();
    }

    @Test
    public void readNext() throws Exception {

    }

    @Test
    public void close() throws Exception {

    }

    @Test
    public void getRecordsLoadedSoFar() throws Exception {

    }

    @Test
    public void getTotalRecords() throws Exception {
        assertTrue("No records found", reader.getRecordsLoadedSoFar() > 0);
    }

}