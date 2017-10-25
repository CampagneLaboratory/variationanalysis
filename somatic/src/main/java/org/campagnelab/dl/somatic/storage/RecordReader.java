package org.campagnelab.dl.somatic.storage;

import org.apache.commons.io.IOUtils;
import org.campagnelab.dl.varanalysis.protobuf.BaseInformationRecords;
import org.campagnelab.goby.baseinfo.SequenceBaseInformationReader;
import org.campagnelab.goby.compression.ChunkCodec;
import org.campagnelab.goby.compression.MessageChunksReader;
import org.campagnelab.goby.exception.GobyRuntimeException;

import java.io.Closeable;
import java.io.File;
import java.io.IOException;
import java.util.Iterator;
import java.util.Properties;
import java.util.Spliterator;
import java.util.function.Consumer;

/**
 * A reader for base information records stored in protobuf format.
 *
 * @author manuele simi
 */
public class RecordReader implements Closeable, RecordIterable<BaseInformationRecords.BaseInformation>,
        RecordReaderI<BaseInformationRecords.BaseInformation> {

    private final SequenceBaseInformationReader reader;

    public RecordReader(String filepath) throws IOException {
        reader = new SequenceBaseInformationReader(filepath);
    }

    /**
     * Creates a reader for the specified range.
     * @param filepath
     * @param startFrom offset in the file to start from (bytes)
     * @param upTo last bytes to read
     * @throws IOException
     */
    public RecordReader(String filepath, long startFrom, long upTo) throws IOException {
        reader = new SequenceBaseInformationReader(startFrom,upTo,filepath);

    }


    /**
     * Reads the next record, if available.
     *
     * @return the record
     * @throws IOException
     */
    public BaseInformationRecords.BaseInformation nextRecord()  {
        try {
            if (reader.hasNext()) return reader.next();
            else return null;
        } catch (GobyRuntimeException e) {
            return null;
        }
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
     * Gets an estimated size of a chunk.
     * @return the size, in bytes
     */
    public long getChunkSize() {
        long fileSize = new File(this.reader.getSourceSbiPath()).length();
        return Math.round( (fileSize / this.getTotalRecords()) * 10000); //10000 is the default number of records per chunk         
    }

    /**
     * Gets an estimated size of a record.
     * @return the size, in bytes
     */
    public long getRecordSize() {
        long fileSize = new File(this.reader.getSourceSbiPath()).length();
        return Math.round(fileSize / this.getTotalRecords()); //10000 is the default number of records per chunk
    }


    /**
     * Gets the current position of the reader.
     * @return the offset of the position, in bytes.
     */
    public long getCurrentReadPosition() {
        return this.reader.getCurrentReadPosition();
    }

        @Override
    public Iterator<BaseInformationRecords.BaseInformation> iterator() {
        return reader.iterator();
    }

    /**
     * Creates a {@link Spliterator} over the elements described by this
     * {@code Iterable}.
     *
     * @return a {@code Spliterator} over the elements described by this
     * {@code Iterable}.
     * @implSpec The default implementation creates an
     * <em><a href="Spliterator.html#binding">early-binding</a></em>
     * spliterator from the iterable's {@code Iterator}.  The spliterator
     * inherits the <em>fail-fast</em> properties of the iterable's iterator.
     * @implNote The default implementation should usually be overridden.  The
     * spliterator returned by the default implementation has poor splitting
     * capabilities, is unsized, and does not report any spliterator
     * characteristics. Implementing classes can nearly always provide a
     * better implementation.
     * @since 1.8
     */
    @Override
    public Spliterator<BaseInformationRecords.BaseInformation> spliterator() {
        long length = new File(this.reader.getSourceSbiPath()).length();
        return new RecordSpliterator(reader.getSourceSbiPath(), 1, length, length);
    }

    public Properties getProperties() {
        return reader.getProperties();
    }


    @Override
    public long numRecords() {
        return getTotalRecords();
    }


    public void readUpTo(long endIndex) {
        this.reader.resetEndOffset(endIndex);
    }
}
