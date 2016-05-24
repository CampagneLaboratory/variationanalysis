package org.campagnelab.dl.varanalysis.storage;


import com.google.protobuf.Message;
import org.apache.commons.io.IOUtils;
import org.apache.hadoop.fs.Path;
import org.apache.parquet.hadoop.ParquetWriter;
import org.apache.parquet.hadoop.metadata.CompressionCodecName;
import org.apache.parquet.proto.ProtoParquetReader;
import org.apache.parquet.proto.ProtoParquetWriter;
import org.campagnelab.dl.varanalysis.protobuf.BaseInformationRecords;

import java.io.Closeable;
import java.io.IOException;
import java.util.List;

/**
 * A writer for base information records in protobuf format.
 *
 * @author manuele simi
 */
public class RecordWriter implements Closeable {

    private final ProtoParquetWriter<BaseInformationRecords.BaseInformation> parquetWriter;
    private static CompressionCodecName compressionCodecName = CompressionCodecName.UNCOMPRESSED;


    public RecordWriter(String file, int blockSize, int pageSize) throws IOException {
       this.parquetWriter =
                new ProtoParquetWriter(new Path(file), BaseInformationRecords.BaseInformation.class, compressionCodecName, blockSize, pageSize);
    }

    public void writeRecord(BaseInformationRecords.BaseInformation record) throws IOException {
        this.parquetWriter.write(record);
    }

    public void writeRecord(List<SampleCountInfo> allCounts, String mutatedBase,
                            String refBase, int refIndex, int position, int index, boolean mutated) throws IOException {
        BaseInformationRecords.BaseInformation.Builder builder = BaseInformationRecords.BaseInformation.newBuilder();
        builder.setMutated(mutated);
        builder.setIndexOfMutatedBase(index);
        builder.setPosition(position);
        builder.setReferenceBase(refBase);
        builder.setReferenceIndex(refIndex);
        builder.setMutatedBase(mutatedBase);
        for (SampleCountInfo counts: allCounts) {
            builder.addCounts(counts.toProtobuf());
        }
        this.parquetWriter.write(builder.build());
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

    public static class SampleCountInfo {

        public boolean matchesReference;

        public String fromSequence, toSequence;

        public int genotypeCountForwardStrand, genotypeCountReverseStrand;


        protected BaseInformationRecords.SampleCountInfo toProtobuf() {
            BaseInformationRecords.SampleCountInfo.Builder builder = BaseInformationRecords.SampleCountInfo.newBuilder();
            builder.setFromSequence(fromSequence);
            builder.setToSequence(toSequence);
            builder.setMatchesReference(matchesReference);
            builder.setGenotypeCountForwardStrand(genotypeCountForwardStrand);
            builder.setGenotypeCountReverseStrand(genotypeCountReverseStrand);
            return builder.build();
        }

    }
}
