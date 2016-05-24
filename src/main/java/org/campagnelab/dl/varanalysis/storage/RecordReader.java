package org.campagnelab.dl.varanalysis.storage;

import org.apache.commons.io.IOUtils;
import org.apache.hadoop.conf.Configuration;
import org.apache.hadoop.fs.Path;
import org.apache.parquet.hadoop.ParquetReader;
import org.apache.parquet.proto.ProtoParquetReader;
import org.campagnelab.dl.varanalysis.protobuf.BaseInformationRecords;

import java.io.*;

/**
 * A reader for base information records stored in protobuf format.
 *
 * @author manuele simi
 */
public class RecordReader implements Closeable {

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
}
