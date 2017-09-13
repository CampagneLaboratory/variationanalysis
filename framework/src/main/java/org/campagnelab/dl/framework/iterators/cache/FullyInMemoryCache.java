package org.campagnelab.dl.framework.iterators.cache;

import it.unimi.dsi.fastutil.objects.ObjectArrayList;
import org.nd4j.linalg.api.ndarray.INDArray;
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
//    private int numDevices;
//    private MagicQueue magicQueue;

    public FullyInMemoryCache(MultiDataSetIterator source) {
        this.source = source;
        if (source.resetSupported()) {
            source.reset();
        }
        reset();
//        if (numDevices >= 2) {
//            magicQueue = new MagicQueue.Builder().setNumberOfBuckets(-1).build();
//        }
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
    public MultiDataSetPreProcessor getPreProcessor() {
        return source.getPreProcessor();
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
    public synchronized void reset() {
//        numDevices = Nd4j.getAffinityManager().getNumberOfDevices();

        index = -1;
        // force traversal and caching of iterator if reset is called before a full traversal:
        if (!sourceIsComplete) {
            cache.clear();
            if (source.resetSupported()) {
                source.reset();
            }
            while (hasNext()) {
                next();
            }
            sourceIsComplete = true;
            reset();

        }
        //  System.out.println("Memory cache size="+cache.size());


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

    int i = 0;

    @Override
    public MultiDataSet next() {
        if (hasNext()) {
            index++;
            if (sourceIsComplete) {
                return cache.get(index);
            } else {

                MultiDataSet multiDataSet = null;
                multiDataSet = source.next();
//                if (numDevices >= 2) {
//                    // try to move the MDS to a random device, so that we won't run out of memory on the first one
//                    // as would happen if all MDS went there (they were built on a cpu)
//
//                    moveToDecide(multiDataSet, (int) (i++ % numDevices));
//                }

                cache.add(multiDataSet);
                return multiDataSet;
            }
        } else {
            throw new NoSuchElementException();
        }
    }

    public static void moveToDecide(MultiDataSet multiDataSet, int deviceIndex) {

        moveArray(deviceIndex, multiDataSet.getFeatures());
        moveArray(deviceIndex, multiDataSet.getLabels());
        moveArray(deviceIndex, multiDataSet.getLabelsMaskArrays());
        moveArray(deviceIndex, multiDataSet.getLabelsMaskArrays());
    }

    public static void moveArray(int deviceIndex, INDArray[] features) {
        if (features == null) return;
        for (INDArray array : features) {
            if (array != null) {
                Nd4j.getAffinityManager().replicateToDevice(deviceIndex, array);
            }
        }
    }

}
