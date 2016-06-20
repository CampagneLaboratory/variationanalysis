package org.campagnelab.dl.varanalysis.storage;


import org.apache.commons.io.FilenameUtils;
import org.apache.commons.io.IOUtils;
import org.apache.hadoop.fs.Path;
import org.apache.parquet.hadoop.metadata.CompressionCodecName;
import org.apache.parquet.proto.ProtoParquetWriter;
import org.campagnelab.dl.varanalysis.protobuf.BaseInformationRecords;

import java.io.Closeable;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.nio.ByteBuffer;
import java.nio.file.FileStore;
import java.nio.file.Files;
//import java.nio.file.Path;

/**
 * A writer for base information records in protobuf format.
 *
 * @author manuele simi
 */
public class RecordWriter implements Closeable {

    private final ProtoParquetWriter<BaseInformationRecords.BaseInformation> parquetWriter;
    private long numRecordsWritten = 0;
    private String path;

    public RecordWriter(java.nio.file.Path path, ProtoParquetWriter<BaseInformationRecords.BaseInformation> parquetWriter) throws IOException {
        this.path=addParqExtension(path.toString());
        this.parquetWriter = parquetWriter;
    }

    public RecordWriter(String file, int blockSize, int pageSize) throws IOException {
        this(new File(addParqExtension(file)).toPath(), new ProtoParquetWriter<BaseInformationRecords.BaseInformation>(new Path(addParqExtension(file)), BaseInformationRecords.BaseInformation.class, CompressionCodecName.UNCOMPRESSED, blockSize, pageSize));

    }

    public RecordWriter(String file, int blockSize, int pageSize, boolean compress) throws IOException {
        this(new File(addParqExtension(file)).toPath(),

                new ProtoParquetWriter<BaseInformationRecords.BaseInformation>(new Path(addParqExtension(file)), BaseInformationRecords.BaseInformation.class,
                        compress ? CompressionCodecName.SNAPPY : CompressionCodecName.UNCOMPRESSED, blockSize, pageSize));
    }

    public void writeRecord(BaseInformationRecords.BaseInformation record) throws IOException {
        this.parquetWriter.write(record);
        numRecordsWritten += 1;
    }

    public static String addParqExtension(String path){
        return  FilenameUtils.removeExtension(path) + ".parquet";
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

        String infoOut = FilenameUtils.removeExtension(path) + ".info";
        new File(infoOut).delete();
        FileWriter infoWriter = new FileWriter(infoOut);
        infoWriter.write(Integer.toString((int)numRecordsWritten));
        infoWriter.close();
        IOUtils.closeQuietly(parquetWriter);
    }

}
