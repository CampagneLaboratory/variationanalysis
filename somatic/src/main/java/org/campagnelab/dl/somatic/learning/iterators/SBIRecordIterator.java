package org.campagnelab.dl.somatic.learning.iterators;

import org.apache.commons.compress.utils.IOUtils;
import org.campagnelab.dl.framework.domains.DomainDescriptor;
import org.campagnelab.dl.varanalysis.protobuf.BaseInformationRecords;
import org.campagnelab.goby.baseinfo.SequenceBaseInformationReader;

import java.io.IOException;

/**
 * An iterator specialized for SBIRecordIterator
 */
public class SBIRecordIterator extends MultiDataSetRecordIterator<BaseInformationRecords.BaseInformationOrBuilder> {
    private String basename;
    private SequenceBaseInformationReader reader;

    public SBIRecordIterator(String inputFilename, int batchSize, DomainDescriptor domainDescriptor) throws IOException {
        super(inputFilename, batchSize, domainDescriptor);
        this.basename = SequenceBaseInformationReader.getBasename(inputFilename);
        reader = new SequenceBaseInformationReader(inputFilename);
    }

    @Override
    public String getBasename() {
        return basename;
    }

    @Override
    long remainingExamples() {
        return reader.getTotalRecords() - reader.getRecordsLoadedSoFar();
    }

    @Override
    public void reset() {
        IOUtils.closeQuietly(reader);
        try {
            reader = new SequenceBaseInformationReader(basename);
        } catch (IOException e) {
            throw new RuntimeException(e);
        }
    }

    @Override
    public boolean hasNext() {
        return hasNextRecord();
    }


    @Override
    public boolean hasNextRecord() {
        return reader.hasNext();
    }

    @Override
    public BaseInformationRecords.BaseInformationOrBuilder nextRecord() {
        return reader.next();
    }

    @Override
    public long getNumRecords() {
        return reader.getTotalRecords();
    }




}
