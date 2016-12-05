package org.campagnelab.dl.framework.iterators.cache;

import org.nd4j.linalg.dataset.api.MultiDataSet;
import org.nd4j.linalg.dataset.api.MultiDataSetPreProcessor;
import org.nd4j.linalg.dataset.api.iterator.MultiDataSetIterator;

import java.util.ArrayList;
import java.util.LinkedList;
import java.util.NoSuchElementException;

/**
 * Fully caches a multi-dataset iterator in memory.
 */
public class FullyInMemoryCache implements MultiDataSetIterator {
    private MultiDataSetIterator source;
    private ArrayList<MultiDataSet> cache = new ArrayList<>();
    private int index = -1;

    private boolean sourceIsComplete;

    public FullyInMemoryCache(MultiDataSetIterator source) {
        this.source = source;
        if (resetSupported()) {
            source.reset();
        }
        reset();

    }

    @Override
    public MultiDataSet next(int num) {
        return next();
    }

    @Override
    public void setPreProcessor(MultiDataSetPreProcessor preProcessor) {
        source.setPreProcessor(preProcessor);
    }

    @Override
    public boolean resetSupported() {
        return true;
    }

    @Override
    public boolean asyncSupported() {
        return false;
    }

    @Override
    public void reset() {
        index = -1;
        // force traversal and caching of iterator if reset is called before a full traversal:
        if (!sourceIsComplete) {
            while (hasNext()) {
                next();
            }
            reset();
        }

    }

    @Override
    public boolean hasNext() {

        if (sourceIsComplete) {
            return index + 1 < cache.size();
        } else {

            final boolean sourceHasNext = source.hasNext();
            if (!sourceHasNext) {
                sourceIsComplete = true;
            }
            return sourceHasNext;
        }
    }

    @Override
    public MultiDataSet next() {
        if (hasNext()) {
            index++;
            if (sourceIsComplete) {
                return cache.get(index);
            } else {

                MultiDataSet o = null;
                o = source.next();
                cache.add(o);
                return o;
            }
        } else {
            throw new NoSuchElementException();
        }
    }
}
