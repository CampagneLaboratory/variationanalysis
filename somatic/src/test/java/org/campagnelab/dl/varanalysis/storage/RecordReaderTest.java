package org.campagnelab.dl.varanalysis.storage;

import org.campagnelab.dl.varanalysis.protobuf.BaseInformationRecords;
import org.junit.After;
import org.junit.Before;
import org.junit.Test;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotNull;

/**
 * Created by mas2182 on 5/24/16.
 * Temporarily disabled for creation of new test data
 */
public class RecordReaderTest {
    private String  filename = "test-data/reader/c1_genotypes_test_proto_VN";

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
    public void readRecords() throws Exception {
        BaseInformationRecords.BaseInformationOrBuilder record = this.reader.nextRecord();
        while (record != null) {
            record = this.reader.nextRecord();
        }
        assertEquals("Records read", 12185, reader.getRecordsLoadedSoFar() );
    }

    @Test
    public void readRecordsWithIterator() throws Exception {
        int numRecordsRead = 0;
        for (BaseInformationRecords.BaseInformationOrBuilder record: this.reader) {
            assertNotNull(record);
            numRecordsRead++;
        }
        assertEquals("Records read", 12185, numRecordsRead);
    }

    @Test
    public void close() throws Exception {

    }

    @Test
    public void getRecordsLoadedSoFar() throws Exception {

    }

    @Test
    public void getTotalRecords() throws Exception {
        assertEquals("Expected records", 12185, reader.getTotalRecords() );
    }

}