package org.campagnelab.dl.framework.training;

import org.campagnelab.dl.framework.iterators.AttachMultiDataSetIterator;
import org.deeplearning4j.datasets.iterator.AsyncMultiDataSetIterator;
import org.nd4j.linalg.dataset.api.iterator.MultiDataSetIterator;

/**
 * A class to wrap a multi-dataset iterator and attach it asynchroneously.
 */
public class WrapInAsyncAttach {
    public static MultiDataSetIterator wrap(MultiDataSetIterator iterator) {
        String prefetchBufferString = System.getProperty("framework.parallelWrapper.prefetchBuffer");
        int prefetchBuffer = prefetchBufferString != null ? Integer.parseInt(prefetchBufferString) : 12;

        //wrap in an async iterator to speed up loading of minibatches to keep the GPU utilized:
        iterator = new AsyncMultiDataSetIterator(iterator, prefetchBuffer);
        // put an attach iterator after the async to attach data to the GPU for all queued mds

        iterator=new AttachMultiDataSetIterator(iterator);
        return iterator;
    }
}
