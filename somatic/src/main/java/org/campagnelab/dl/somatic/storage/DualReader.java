package org.campagnelab.dl.somatic.storage;

import org.apache.commons.compress.utils.IOUtils;
import org.campagnelab.dl.varanalysis.protobuf.BaseInformationRecords;

import java.io.Closeable;
import java.io.IOException;
import java.util.NoSuchElementException;

/**
 * Helper class to iterate over two readers that are aligned by genomic position.
 * Created by fac2003 on 5/26/16.
 *
 * @author Fabien Campagne
 */
public class DualReader implements Closeable {
    private final long totalRecords;
    private RecordReader readerA;
    private RecordReader readerB;
    private BaseInformationRecords.BaseInformation recordA;
    private BaseInformationRecords.BaseInformation recordB;
    private boolean loaded;

    public DualReader(String filename1, String filename2) throws IOException {
        readerA = new RecordReader(filename1);
        readerB = new RecordReader(filename2);
        totalRecords=readerA.getTotalRecords();
        assert readerA.getTotalRecords() == readerB.getTotalRecords() : " readers must have the same number of records";
    }

    public long getTotalRecords() {
        return totalRecords ;
    }

    public boolean hasNext() throws IOException {

        if (loaded) return true;
        else {
            recordA = readerA.nextRecord();
            recordB = readerB.nextRecord();
            if (recordA==null || recordB==null) {
                loaded=false;
                return false;
            }
            assert recordA.getReferenceIndex() == recordB.getReferenceIndex() : "reference indices must match between readers.";
            assert recordA.getPosition() == recordB.getPosition() : "position must match between readers.";
        }
        loaded = recordA != null && recordB != null;
        return loaded;
    }

    /**
     * Return the next element from the first reader.
     *
     * @return
     */
    public BaseInformationRecords.BaseInformation first() {
        if (!loaded) throw new NoSuchElementException();
        else {
            return recordA;
        }
    }

    /**
     * Return the next element from the second reader.
     *
     * @return
     */
    public BaseInformationRecords.BaseInformation second() {
        if (!loaded) throw new NoSuchElementException();
        else {
            loaded=false;
            return recordB;

        }
    }

    public void close() {
        IOUtils.closeQuietly(readerA);
        IOUtils.closeQuietly(readerB);

    }
}
