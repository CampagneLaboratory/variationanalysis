package org.campagnelab.dl.varanalysis.learning;

import java.io.Closeable;

/**
 * Interface for readers over training records.
 */
public interface RecordReaderI<RecordType> extends Closeable, Iterable<RecordType> {
    /**
     * The number of records that can be read with this reader.
     * @return number of records
     */
    long numRecords();

}
