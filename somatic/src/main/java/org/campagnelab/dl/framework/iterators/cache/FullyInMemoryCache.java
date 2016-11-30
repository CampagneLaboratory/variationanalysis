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
    private int cachedIndex = -1;
    private int index = -1;

    public FullyInMemoryCache(MultiDataSetIterator source) {
        this.source = source;

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

    }

    @Override
    public boolean hasNext() {

        if (index < cachedIndex) {
            return true;
        } else {
            return source.hasNext();
        }
    }

    @Override
    public MultiDataSet next() {
        if (hasNext()) {
            index++;
            if (index >= cachedIndex) {

                MultiDataSet o = source.next();
                cache.add(o);
                cachedIndex++;
            }
            return cache.get(index);
        } else {
            throw new NoSuchElementException();
        }
    }
}
