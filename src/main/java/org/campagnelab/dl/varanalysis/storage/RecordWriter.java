package org.campagnelab.dl.varanalysis.storage;


import org.apache.commons.io.IOUtils;
import org.apache.hadoop.fs.Path;
import org.apache.parquet.hadoop.metadata.CompressionCodecName;
import org.apache.parquet.proto.ProtoParquetWriter;
import org.campagnelab.dl.varanalysis.protobuf.BaseInformationRecords;

import java.io.Closeable;
import java.io.File;
import java.io.IOException;
import java.nio.ByteBuffer;
import java.nio.file.FileStore;
import java.nio.file.Files;
//import java.nio.file.Path;
import java.nio.file.attribute.UserDefinedFileAttributeView;

/**
 * A writer for base information records in protobuf format.
 *
 * @author manuele simi
 */
public class RecordWriter implements Closeable {

    private final ProtoParquetWriter<BaseInformationRecords.BaseInformation> parquetWriter;
    private long numRecordsWritten = 0;
    private FileStore store;
    private UserDefinedFileAttributeView view;

    public RecordWriter(java.nio.file.Path path, ProtoParquetWriter<BaseInformationRecords.BaseInformation> parquetWriter) throws IOException {
        this.parquetWriter = parquetWriter;
        // check that user defined attributes are supported by the file store
        store = Files.getFileStore(path);
        if (store.supportsFileAttributeView(UserDefinedFileAttributeView.class)) {
            view = Files.
                    getFileAttributeView(path, UserDefinedFileAttributeView.class);
            System.err.println("File system supports attributes: YES");
        }else {
            System.err.println("File system supports attributes: NO");
        }


    }

    public RecordWriter(String file, int blockSize, int pageSize) throws IOException {

        this(new File(file).toPath(), new ProtoParquetWriter<BaseInformationRecords.BaseInformation>(new Path(file), BaseInformationRecords.BaseInformation.class, CompressionCodecName.UNCOMPRESSED, blockSize, pageSize));

    }

    public RecordWriter(String file, int blockSize, int pageSize, boolean compress) throws IOException {
        this(new File(file).toPath(),

                new ProtoParquetWriter<BaseInformationRecords.BaseInformation>(new Path(file), BaseInformationRecords.BaseInformation.class,
                        compress ? CompressionCodecName.SNAPPY : CompressionCodecName.UNCOMPRESSED, blockSize, pageSize));

    }

    public void writeRecord(BaseInformationRecords.BaseInformation record) throws IOException {
        this.parquetWriter.write(record);
        numRecordsWritten += 1;
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
        if (view!=null) {
            //file system supports user defined file attributes.
            // we store the total number of records written to the file, so that we avoid scanning this file just to
            // determine this number upon read:
            ByteBuffer b=ByteBuffer.allocate(8);
            b.asLongBuffer().put(numRecordsWritten);

            view.write("numRecordsWritten", b);
        }
        IOUtils.closeQuietly(parquetWriter);
    }

}
