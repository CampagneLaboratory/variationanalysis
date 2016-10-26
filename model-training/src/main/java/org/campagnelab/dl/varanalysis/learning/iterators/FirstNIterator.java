package org.campagnelab.dl.varanalysis.learning.iterators;

import org.nd4j.linalg.dataset.DataSet;
import org.nd4j.linalg.dataset.api.DataSetPreProcessor;
import org.nd4j.linalg.dataset.api.iterator.DataSetIterator;

import java.util.List;
import java.util.NoSuchElementException;

/**
 * Created by fac2003 on 7/14/16.
 */
public class FirstNIterator implements DataSetIterator {
    private final int N;
    DataSetIterator delegate;
    int index;

    @Override
    public DataSet next(int num) {
        if (index >= N) {
            throw new NoSuchElementException();
        }
        final DataSet next = delegate.next(num);
        index += next.numExamples();
        return next;
    }

    @Override
    public int totalExamples() {
        return Math.min(N, delegate.totalExamples());
    }

    @Override
    public int inputColumns() {
        return delegate.inputColumns();
    }

    @Override
    public int totalOutcomes() {
        return Math.min(N, delegate.totalOutcomes());
    }

    @Override
    public boolean resetSupported() {
        return true;
    }


    public boolean asyncSupported() {
        return false;
    }

    @Override
    public void reset() {
        index = 0;
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

        return Math.min(N, delegate.numExamples());
    }

    @Override
    public void setPreProcessor(DataSetPreProcessor preProcessor) {
        delegate.setPreProcessor(preProcessor);
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
        if (index >= N) {
            return false;
        }
        return delegate.hasNext();
    }

    @Override
    public DataSet next() {
        if (index >= N) {
            throw new NoSuchElementException();
        }

        DataSet next = delegate.next();
        index += next.numExamples();
        return next;
    }

    public FirstNIterator(DataSetIterator delegate, int n) {
        this.delegate = delegate;
        this.N = n;
    }
}
