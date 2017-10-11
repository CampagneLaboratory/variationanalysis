package org.campagnelab.dl.somatic.storage;

import org.campagnelab.dl.varanalysis.protobuf.SegmentInformationRecords;
import org.campagnelab.goby.baseinfo.SequenceSegmentInformationWriter;

import java.io.Closeable;
import java.io.IOException;

/**
 * A writer for sequence segment information records in protobuf format.
 *
 * @author manuele simi
 */
public class SegmentWriter implements Closeable {

    private SequenceSegmentInformationWriter writer;

    public SegmentWriter(String file, int numEntriesPerChunk) throws IOException {
        writer = new SequenceSegmentInformationWriter(file);
        writer.setNumEntriesPerChunk(numEntriesPerChunk);
    }

    public SegmentWriter(String file) throws IOException {
        writer = new SequenceSegmentInformationWriter(file);
    }

    public void writeRecord(SegmentInformationRecords.SegmentInformation record) throws IOException {
        writer.appendEntry(record);
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
