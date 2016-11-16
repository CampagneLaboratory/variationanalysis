package org.campagnelab.dl.somatic.learning.iterators;

import org.nd4j.linalg.dataset.DataSet;
import org.nd4j.linalg.dataset.api.DataSetPreProcessor;
import org.nd4j.linalg.dataset.api.iterator.DataSetIterator;

import java.util.List;

/**
 * An adapter to provide a basename for a CachingDatasetIterator.
 * Created by fac2003 on 11/2/16.
 */
public class NamedCachingDataSetIterator implements NamedDataSetIterator {
    @Override
    public DataSet next(int i) {
        return delegate.next(i);
    }

    @Override
    public int totalExamples() {
        return delegate.totalExamples();
    }

    @Override
    public int inputColumns() {
        return delegate.inputColumns();
    }

    @Override
    public int totalOutcomes() {
        return delegate.totalOutcomes();
    }

    @Override
    public boolean resetSupported() {
        return delegate.resetSupported();
    }


    public boolean asyncSupported() {
        return false;
    }

    @Override
    public void reset() {
        delegate.reset();
    }

    @Override
    public int batch() {
        return delegate.batch();
    }

    @Override
    public int cursor() {
        return delegate.cursor();
    }

    @Override
    public int numExamples() {
        return delegate.numExamples();
    }

    @Override
    public void setPreProcessor(DataSetPreProcessor dataSetPreProcessor) {
        delegate.setPreProcessor(dataSetPreProcessor);
    }

    @Override
    public DataSetPreProcessor getPreProcessor() {
        return delegate.getPreProcessor();
    }

    @Override
    public List<String> getLabels() {
        return delegate.getLabels();
    }

    @Override
    public boolean hasNext() {
        return delegate.hasNext();
    }

    @Override
    public DataSet next() {
        return delegate.next();
    }

    private DataSetIterator delegate;
    private String basename;

    public NamedCachingDataSetIterator(DataSetIterator delegate, String basename) {
        this.delegate = delegate;
    }

    @Override
    public String getBasename() {
        return basename;
    }
}
