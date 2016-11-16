package org.campagnelab.dl.somatic.storage;


import org.apache.commons.io.FilenameUtils;
import org.campagnelab.dl.varanalysis.protobuf.BaseInformationRecords;
import org.campagnelab.goby.baseinfo.SequenceBaseInformationWriter;

import java.io.Closeable;
import java.io.IOException;

/**
 * A writer for base information records in protobuf format.
 *
 * @author manuele simi
 */
public class RecordWriter implements Closeable {

   private SequenceBaseInformationWriter writer;
    public RecordWriter(String file, int numEntriesPerChunk) throws IOException {
        writer = new SequenceBaseInformationWriter(file);
        writer.setNumEntriesPerChunk(numEntriesPerChunk);
    }

    public RecordWriter(String file) throws IOException {
        writer = new SequenceBaseInformationWriter(file);
    }

    public void writeRecord(BaseInformationRecords.BaseInformation record) throws IOException {
        writer.appendEntry(record);
    }

    public static String addParqExtension(String path) {
        return FilenameUtils.removeExtension(path) + ".parquet";
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

        writer.close();
    }

}
