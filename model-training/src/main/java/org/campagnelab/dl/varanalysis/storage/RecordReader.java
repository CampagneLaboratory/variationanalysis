package org.campagnelab.dl.varanalysis.storage;

import com.codahale.metrics.MetricRegistryListener;
import com.google.protobuf.CodedInputStream;
import it.unimi.dsi.fastutil.ints.Int2IntArrayMap;
import it.unimi.dsi.fastutil.ints.Int2IntMap;
import it.unimi.dsi.fastutil.ints.IntArrayList;
import it.unimi.dsi.fastutil.objects.ObjectArrayList;
import it.unimi.dsi.logging.ProgressLogger;
import org.apache.commons.io.FilenameUtils;
import org.apache.commons.io.IOUtils;
import org.apache.hadoop.conf.Configuration;
import org.apache.hadoop.fs.Path;
import org.apache.parquet.hadoop.ParquetReader;
import org.apache.parquet.hadoop.api.ReadSupport;
import org.apache.parquet.io.api.RecordMaterializer;
import org.apache.parquet.proto.ProtoParquetReader;
import org.apache.parquet.schema.MessageType;
import org.campagnelab.dl.varanalysis.protobuf.BaseInformationRecords;

import java.io.*;
import java.nio.ByteBuffer;
import java.nio.file.FileStore;
import java.nio.file.Files;
import java.nio.file.attribute.UserDefinedFileAttributeView;
import java.util.*;
import java.util.function.Consumer;

/**
 * A reader for base information records stored in protobuf format.
 *
 * @author manuele simi
 */
public class RecordReader implements Closeable, RecordIterable {

    private final String filepath;
    private final ParquetReader<BaseInformationRecords.BaseInformation> pqReader;
    private UserDefinedFileAttributeView view;
    private long recordLoadedSoFar = 0;
    private long totalRecords = -1;

    public RecordReader(String filepath) throws IOException {
        this.filepath=RecordWriter.addParqExtension(filepath);
        this.countRecords();

        ProtoParquetReader.Builder<BaseInformationRecords.BaseInformation> pqBuilder = ProtoParquetReader.builder(new Path(filepath));
        Configuration conf = new Configuration();
        conf.set("parquet.proto.class", BaseInformationRecords.BaseInformation.class.getCanonicalName());
        pqBuilder.withConf(conf);
        pqReader = pqBuilder.build();
    }

    private void countRecords() throws IOException {
        String countPath = FilenameUtils.removeExtension(filepath) + ".info";
        File countFile = new File(countPath);
        if (countFile.exists()) {
            BufferedReader countReader = new BufferedReader(new FileReader(countFile));
            totalRecords = Integer.valueOf(countReader.readLine());
        } else {
            System.out.println("count file not found: " + filepath.substring(0, filepath.length() - 8) + ".info");
            totalRecords = 0;
            ProtoParquetReader.Builder<BaseInformationRecords.BaseInformation> pqBuilder = ProtoParquetReader.builder(new Path(this.filepath));
            Configuration conf = new Configuration();
            conf.set("parquet.proto.class", BaseInformationRecords.BaseInformation.class.getCanonicalName());
            pqBuilder.withConf(conf);
            ParquetReader<BaseInformationRecords.BaseInformation> localReader = pqBuilder.build();
            while (localReader.read() != null) this.totalRecords++;
            IOUtils.closeQuietly(localReader);
        }

    }

    /**
     * Reads the next record, if available.
     *
     * @return the record
     * @throws IOException
     */
    public BaseInformationRecords.BaseInformation nextRecord() throws IOException {
        BaseInformationRecords.BaseInformationOrBuilder record = pqReader.read();

        if (record == null) {
            return null;
        }
        recordLoadedSoFar++;

        if (record instanceof BaseInformationRecords.BaseInformation.Builder) {
            return ((BaseInformationRecords.BaseInformation.Builder) record).build();
        } else {
            assert false : "Cannot build record.";
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
    public void forEach(Consumer<? super BaseInformationRecords.BaseInformation> action) {

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
        return null;
    }


}
