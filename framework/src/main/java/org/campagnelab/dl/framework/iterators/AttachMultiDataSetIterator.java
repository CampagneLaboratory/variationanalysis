package org.campagnelab.dl.framework.iterators;

import org.nd4j.linalg.dataset.api.MultiDataSet;
import org.nd4j.linalg.dataset.api.MultiDataSetPreProcessor;
import org.nd4j.linalg.dataset.api.iterator.MultiDataSetIterator;

/**
 * An iterator that attaches multidatasets before returning them in next.
 * Created by fac2003 on 11/4/17.
 */
public class AttachMultiDataSetIterator implements MultiDataSetIterator {
    MultiDataSetIterator delegate;
    MultiDataSetPreProcessor preProcessor;

    @Override
    public MultiDataSet next(int num) {
        MultiDataSet next = delegate.next(num);
        if (preProcessor != null) preProcessor.preProcess(next);
        MDSHelper.attach(next);
        return next;
    }

    @Override
    /** Set the preprocessor of this iterator. */
    public void setPreProcessor(MultiDataSetPreProcessor preProcessor) {
        this.preProcessor = preProcessor;
    }

    @Override
    public MultiDataSetPreProcessor getPreProcessor() {
        return preProcessor;
    }

    @Override
    public boolean resetSupported() {
        return delegate.resetSupported();
    }

    @Override
    public boolean asyncSupported() {
        return delegate.asyncSupported();
    }

    @Override
    public void reset() {
        delegate.reset();
    }

    @Override
    public boolean hasNext() {
        return delegate.hasNext();
    }

    @Override
    public MultiDataSet next() {
        MultiDataSet next = delegate.next();
        if (preProcessor != null) preProcessor.preProcess(next);
        MDSHelper.attach(next);
        return next;
    }

    public AttachMultiDataSetIterator(MultiDataSetIterator delegate) {
        this.delegate = delegate;
    }
}
