package org.campagnelab.dl.varanalysis.storage;


import org.apache.commons.io.IOUtils;
import org.apache.hadoop.fs.Path;
import org.apache.parquet.hadoop.metadata.CompressionCodecName;
import org.apache.parquet.proto.ProtoParquetWriter;
import org.campagnelab.dl.varanalysis.protobuf.BaseInformationRecords;

import java.io.Closeable;
import java.io.IOException;

/**
 * A writer for base information records in protobuf format.
 *
 * @author manuele simi
 */
public class RecordWriter implements Closeable {

    private final ProtoParquetWriter<BaseInformationRecords.BaseInformation> parquetWriter;

    public RecordWriter(String file, int blockSize, int pageSize) throws IOException {
        this.parquetWriter =
                new ProtoParquetWriter<BaseInformationRecords.BaseInformation>(new Path(file), BaseInformationRecords.BaseInformation.class, CompressionCodecName.UNCOMPRESSED, blockSize, pageSize);
    }

    public RecordWriter(String file, int blockSize, int pageSize, boolean compress) throws IOException {
        this.parquetWriter =
                new ProtoParquetWriter<BaseInformationRecords.BaseInformation>(new Path(file), BaseInformationRecords.BaseInformation.class,
                        compress ? CompressionCodecName.SNAPPY : CompressionCodecName.UNCOMPRESSED, blockSize, pageSize);
    }

    public void writeRecord(BaseInformationRecords.BaseInformation record) throws IOException {
        this.parquetWriter.write(record);
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
        IOUtils.closeQuietly(parquetWriter);
    }

}
