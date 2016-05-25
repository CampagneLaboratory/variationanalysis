package org.campagnelab.dl.varanalysis.storage;

import org.apache.commons.io.IOUtils;
import org.apache.hadoop.conf.Configuration;
import org.apache.hadoop.fs.Path;
import org.apache.parquet.hadoop.ParquetReader;
import org.apache.parquet.proto.ProtoParquetReader;
import org.campagnelab.dl.varanalysis.protobuf.BaseInformationRecords;

import java.io.*;
import java.util.Iterator;
import java.util.Spliterator;
import java.util.function.Consumer;

/**
 * A reader for base information records stored in protobuf format.
 *
 * @author manuele simi
 */
public class RecordReader implements Closeable, RecordIterable {

    private final String filepath;
    private final ParquetReader<BaseInformationRecords.BaseInformation> pqReader;
    private long recordLoadedSoFar = 0;
    private long totalRecords = 0;

    public RecordReader(String filepath) throws IOException {
        this.filepath = filepath;
        this.countRecords();
        ProtoParquetReader.Builder<BaseInformationRecords.BaseInformation> pqBuilder = ProtoParquetReader.builder(new Path(filepath));
        Configuration conf = new Configuration();
        conf.set("parquet.proto.class", BaseInformationRecords.BaseInformation.class.getCanonicalName());
        pqBuilder.withConf(conf);
        pqReader = pqBuilder.build();
    }

    private void countRecords() throws IOException {
        ProtoParquetReader.Builder<BaseInformationRecords.BaseInformation> pqBuilder = ProtoParquetReader.builder(new Path(this.filepath));
        Configuration conf = new Configuration();
                conf.set("parquet.proto.class", BaseInformationRecords.BaseInformation.class.getCanonicalName());
        pqBuilder.withConf(conf);
        ParquetReader<BaseInformationRecords.BaseInformation> localReader = pqBuilder.build();
        while (localReader.read() != null) this.totalRecords++;
        IOUtils.closeQuietly(localReader);
    }

    /**
     * Reads the next record, if available.
     *
     * @return the record
     * @throws IOException
     */
    public BaseInformationRecords.BaseInformationOrBuilder nextRecord() throws IOException {
        BaseInformationRecords.BaseInformationOrBuilder record = pqReader.read();
        if (record == null)
            return null;
        recordLoadedSoFar++;
        return record;
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
        IOUtils.closeQuietly(pqReader);
    }

    /**
     * Gets the number of records read so far.
     *
     * @return records loaded
     */
    public long getRecordsLoadedSoFar() {
        return this.recordLoadedSoFar;
    }

    /**
     * Gets the total number of records.
     *
     * @return total records
     */
    public long getTotalRecords() {
        return this.totalRecords;
    }

    /**
     * Returns an iterator over elements of type {@code T}.
     *
     * @return an Iterator.
     */
    @Override
    public RecordIterator iterator() {
        ProtoParquetReader.Builder<BaseInformationRecords.BaseInformation> pqBuilder = ProtoParquetReader.builder(new Path(filepath));
        Configuration conf = new Configuration();
        conf.set("parquet.proto.class", BaseInformationRecords.BaseInformation.class.getCanonicalName());
        pqBuilder.withConf(conf);
        try {
            RecordIterator iterator = new RecordIterator(pqBuilder.build());
            iterator.setTotalRecords(this.totalRecords);
            return iterator;
        } catch (IOException e) {
            return null;
        }
    }

    /**
     * Performs the given action for each element of the {@code Iterable}
     * until all elements have been processed or the action throws an
     * exception.  Unless otherwise specified by the implementing class,
     * actions are performed in the order of iteration (if an iteration order
     * is specified).  Exceptions thrown by the action are relayed to the
     * caller.
     *
     * @param action The action to be performed for each element
     * @throws NullPointerException if the specified action is null
     * @implSpec <p>The default implementation behaves as if:
     * <pre>{@code
     *     for (T t : this)
     *         action.accept(t);
     * }</pre>
     * @since 1.8
     */
    @Override
    public void forEach(Consumer<? super BaseInformationRecords.BaseInformationOrBuilder> action) {

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
    public Spliterator<BaseInformationRecords.BaseInformationOrBuilder> spliterator() {
        return null;
    }
}
