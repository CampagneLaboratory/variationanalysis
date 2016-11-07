package org.campagnelab.dl.varanalysis.learning.iterators;

import org.nd4j.linalg.dataset.api.MultiDataSet;
import org.nd4j.linalg.dataset.api.MultiDataSetPreProcessor;
import org.nd4j.linalg.dataset.api.iterator.MultiDataSetIterator;

/**
 * Created by fac2003 on 11/7/16.
 */
public class MultiDataSetIteratorAdapter implements MultiDataSetIterator {


    @Override
    public MultiDataSet next(int num) {
        return null;
    }

    @Override
    public void setPreProcessor(MultiDataSetPreProcessor preProcessor) {

    }

    @Override
    public boolean resetSupported() {
        return false;
    }

    @Override
    public boolean asyncSupported() {
        return false;
    }

    @Override
    public void reset() {

    }

    @Override
    public boolean hasNext() {
        return false;
    }

    @Override
    public MultiDataSet next() {
        return null;
    }
}
