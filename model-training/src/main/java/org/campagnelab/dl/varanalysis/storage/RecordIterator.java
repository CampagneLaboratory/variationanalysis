package org.campagnelab.dl.varanalysis.storage;

import org.apache.commons.compress.utils.IOUtils;
import org.apache.parquet.hadoop.ParquetReader;
import org.campagnelab.dl.varanalysis.protobuf.BaseInformationRecords;

import java.io.IOException;
import java.util.Iterator;

/**
 * An Iterator over {@link ParquetReader}.
 *
 * @author Manuele Simi
 */
public class RecordIterator implements Iterator<BaseInformationRecords.BaseInformation> {

    private long totalRecords = 0;
    private long recordLoadedSoFar = 0;
    private final ParquetReader<BaseInformationRecords.BaseInformation> reader;

    protected RecordIterator(ParquetReader<BaseInformationRecords.BaseInformation> reader) {
      this.reader = reader;
    }

    /**
     * Expected records to be read by the iterator.
     * @param totalRecords
     */
    protected void setTotalRecords(long totalRecords) {
        this.totalRecords = totalRecords;
    }

    @Override
    public boolean hasNext() {
        return (recordLoadedSoFar < totalRecords);
    }

    @Override
    public BaseInformationRecords.BaseInformation next() {
        BaseInformationRecords.BaseInformationOrBuilder record = null;
        try {
            record = reader.read();
            recordLoadedSoFar++;
        } catch (IOException e) {
        } finally {
            if (record instanceof BaseInformationRecords.BaseInformation.Builder) {
                return  ((BaseInformationRecords.BaseInformation.Builder) record).build();
            }else {
                assert false:"Cannot build record.";
                return null;
            }
        }
    }

    /**
     * Called by the garbage collector on an object when garbage collection
     * determines that there are no more references to the iterator.
     */
    @Override
    protected void finalize() throws Throwable {
        super.finalize();
        IOUtils.closeQuietly(this.reader);
    }
}
