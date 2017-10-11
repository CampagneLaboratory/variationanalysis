package org.campagnelab.dl.somatic.storage;

import org.apache.commons.io.IOUtils;
import org.campagnelab.dl.varanalysis.protobuf.SegmentInformationRecords;
import org.campagnelab.goby.baseinfo.SequenceSegmentInformationReader;
import org.campagnelab.goby.exception.GobyRuntimeException;

import java.io.Closeable;
import java.io.IOException;
import java.util.Iterator;
import java.util.Properties;

/**
 * A reader for segment information records stored in protobuf format.
 *
 * @author manuele
 */
public class SegmentReader implements Closeable, RecordIterable<SegmentInformationRecords.SegmentInformation>,
        RecordReaderI<SegmentInformationRecords.SegmentInformation> {

    private SequenceSegmentInformationReader reader;

    public SegmentReader(String filepath) throws IOException {
        reader = new SequenceSegmentInformationReader(filepath);
    }
    /**
     * The number of records that can be read with this reader.
     *
     * @return number of records
     */
    @Override
    public long numRecords() {
        return getTotalRecords();
    }
    /**
     * Reads the next segment, if available.
     *
     * @return the record
     * @throws IOException
     */
    public SegmentInformationRecords.SegmentInformation nextSegment() throws IOException {
        try {
            if (reader.hasNext()) return reader.next();
            else return null;
        } catch (GobyRuntimeException e) {
            return null;
        }
    }

    public Properties getProperties() {
        return reader.getProperties();
    }

    /**
     * Gets the number of records read so far.
     *
     * @return records loaded
     */
    public long getRecordsLoadedSoFar() {
        return reader.getRecordsLoadedSoFar();
    }

    /**
     * Gets the total number of records.
     *
     * @return total records
     */
    public long getTotalRecords() {
        return reader.getTotalRecords();
    }

    /**
     * Closes this stream and releases any system resources associated
     * with it. If the stream is already closed then invoking this
     * method has no effect.
     * <p>
     * <p> As noted in {@link AutoCloseable#close()}, cases where the
     * close may fail require careful attention. It is strongly advised
     * to relinquish the underlying resources and to internally
     * <em>mark</em> the {@code Closeable} as closed, prior to throwing
     * the {@code IOException}.
     *
     * @throws IOException if an I/O error occurs
     */
    @Override
    public void close() throws IOException {
        IOUtils.closeQuietly(reader);
    }

    /**
     * Returns an iterator over elements of type {@code T}.
     *
     * @return an Iterator.
     */
    @Override
    public Iterator<SegmentInformationRecords.SegmentInformation> iterator() {
        return reader.iterator();
    }
}
