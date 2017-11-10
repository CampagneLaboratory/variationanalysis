package org.campagnelab.dl.framework.iterators;

import org.nd4j.linalg.dataset.api.MultiDataSet;
import org.nd4j.linalg.dataset.api.MultiDataSetPreProcessor;
import org.nd4j.linalg.dataset.api.iterator.MultiDataSetIterator;

/**
 * An iterator that is the same as the source, but returns asyncSupported(). Useful to trick ParallelWrapper
 * when we have already installed an async iterator (which returns asyncSupported=false).
 */
public class ForceAsync implements MultiDataSetIterator {
    MultiDataSetIterator delegate;

    @Override
    public MultiDataSet next(int num) {
        return delegate.next(num);
    }

    @Override
    public void setPreProcessor(MultiDataSetPreProcessor preProcessor) {
        delegate.setPreProcessor(preProcessor);
    }

    @Override
    public MultiDataSetPreProcessor getPreProcessor() {
        return delegate.getPreProcessor();
    }

    @Override
    public boolean resetSupported() {
        return delegate.resetSupported();
    }

    @Override
    public boolean asyncSupported() {
        return true;
    }

    @Override
    public void reset() {
        delegate.reset();
    }

    public ForceAsync(MultiDataSetIterator iterator) {
        this.delegate=iterator;
    }

    @Override
    public boolean hasNext() {
        return delegate.hasNext();
    }

    @Override
    public MultiDataSet next() {
        return delegate.next();
    }
}
