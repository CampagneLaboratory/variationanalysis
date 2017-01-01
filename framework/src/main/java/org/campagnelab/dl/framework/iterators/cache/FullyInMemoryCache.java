package org.campagnelab.dl.framework.iterators.cache;

import it.unimi.dsi.fastutil.objects.ObjectArrayList;
import org.nd4j.linalg.api.concurrency.AffinityManager;
import org.nd4j.linalg.dataset.api.MultiDataSet;
import org.nd4j.linalg.dataset.api.MultiDataSetPreProcessor;
import org.nd4j.linalg.dataset.api.iterator.MultiDataSetIterator;
import org.nd4j.linalg.factory.Nd4j;

import java.util.NoSuchElementException;

/**
 * Fully caches a multi-dataset iterator in memory.
 */
public class FullyInMemoryCache implements MultiDataSetIterator {
    private MultiDataSetIterator source;
    private ObjectArrayList<MultiDataSet> cache = new ObjectArrayList<>();
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
            cache.clear();
            while (hasNext()) {
                next();
            }
            reset();
            System.out.println("Memory cache size="+cache.size());
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

                MultiDataSet multiDataSet = null;
                multiDataSet = source.next();
               // assignToHostMemory(multiDataSet);
                cache.add(multiDataSet);
                return multiDataSet;
            }
        } else {
            throw new NoSuchElementException();
        }
    }

    private void assignToHostMemory(MultiDataSet dataSet) {

        for (int i = 0; i < dataSet.numFeatureArrays(); i++) {
            Nd4j.getAffinityManager().tagLocation(dataSet.getFeatures(i), AffinityManager.Location.HOST);
            Nd4j.getAffinityManager().tagLocation(dataSet.getFeaturesMaskArray(i), AffinityManager.Location.HOST);
        }
        for (int i = 0; i < dataSet.numLabelsArrays(); i++) {
            Nd4j.getAffinityManager().tagLocation(dataSet.getLabels(i), AffinityManager.Location.HOST);
            Nd4j.getAffinityManager().tagLocation(dataSet.getLabelsMaskArray(i), AffinityManager.Location.HOST);
        }
    }
}
