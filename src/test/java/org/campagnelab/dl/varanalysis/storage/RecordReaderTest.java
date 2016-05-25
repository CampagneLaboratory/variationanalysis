package org.campagnelab.dl.varanalysis.storage;

import org.campagnelab.dl.varanalysis.protobuf.BaseInformationRecords;
import org.junit.After;
import org.junit.Before;
import org.junit.Test;

import static org.junit.Assert.*;

/**
 * Created by mas2182 on 5/24/16.
 */
public class RecordReaderTest {
    private String  filename = "test-results/genotypes_mutated_protofbuf.parquet";

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
    public void readrecords() throws Exception {
        BaseInformationRecords.BaseInformationOrBuilder record = this.reader.nextRecord();
        while (record != null) {
            record = this.reader.nextRecord();
        }
        assertEquals("Records read", 100, reader.getRecordsLoadedSoFar() );
    }

    @Test
    public void readrecordsWithIterator() throws Exception {
        int numRecordsRead = 0;
        for (BaseInformationRecords.BaseInformationOrBuilder record: this.reader) {
            assertNotNull(record);
            numRecordsRead++;
        }
        assertEquals("Records read", 100, numRecordsRead);
    }

    @Test
    public void close() throws Exception {

    }

    @Test
    public void getRecordsLoadedSoFar() throws Exception {

    }

    @Test
    public void getTotalRecords() throws Exception {
        assertEquals("Expected records", 100, reader.getTotalRecords() );
    }

}