package org.campagnelab.dl.varanalysis.learning.iterators;

import it.unimi.dsi.fastutil.objects.ObjectArrayList;
import org.campagnelab.dl.model.utils.mappers.FeatureMapper;
import org.campagnelab.dl.model.utils.mappers.LabelMapper;
import org.campagnelab.dl.varanalysis.protobuf.BaseInformationRecords;
import org.deeplearning4j.datasets.iterator.DataSetIterator;
import org.nd4j.linalg.api.ndarray.INDArray;
import org.nd4j.linalg.dataset.DataSet;
import org.nd4j.linalg.dataset.api.DataSetPreProcessor;
import org.nd4j.linalg.factory.Nd4j;

import java.io.IOException;
import java.util.List;
import java.util.NoSuchElementException;

/**
 * Created by rct66 on 6/13/16.
 */
public class BaseInformationConcatIterator extends BaseInformationIterator implements DataSetIterator {

    private List<BaseInformationIterator> baseIters;
    private int readerIndex = 0;

    public BaseInformationConcatIterator(List<BaseInformationIterator> iterators, int batchSize, FeatureMapper featureMapper, LabelMapper labelMapper) throws IOException {
        super(featureMapper,labelMapper);
        this.batchSize = batchSize;
        this.baseIters = iterators;
        for (BaseInformationIterator iter : iterators) {
            this.totalExamples+=iter.totalExamples();
        }

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
        for (BaseInformationIterator iter : baseIters) {
            iter.reset();
        }

        cursor = 0;
        readerIndex = 0;
        nextPosRecord = null;
        //    System.out.println("reset called");
    }


    @Override
    public void setPreProcessor(DataSetPreProcessor dataSetPreProcessor) {
    }

    @Override
    public DataSetPreProcessor getPreProcessor() {
        return null;
    }

    @Override
    public List<String> getLabels() {
        throw new UnsupportedOperationException("Not implemented for the entire dataset.");
    }

    @Override
    public boolean hasNext() {
        if (nextPosRecord != null) {
            return true;
        }
        if (readerIndex >= baseIters.size()){
            return false;
        }
        try {
            this.nextPosRecord = baseIters.get(readerIndex).nextRecord();
        } catch ( NoSuchElementException e){
            readerIndex++;
            return hasNext();

        }
        return nextPosRecord != null;
    }


    @Override
    public DataSet next() {
        if (hasNext()) {
            return next(batchSize);
        } else throw new NoSuchElementException();
    }

    /**
     * Check if there is a next single PosRecord.
     *
     * @return
     */
    public boolean hasNextRecord() {
        if (nextPosRecord != null) {
            return true;
        }
        if (readerIndex >= baseIters.size()){
            return false;
        }
        try {
            this.nextPosRecord = baseIters.get(readerIndex).nextRecord();
        } catch ( NoSuchElementException e){
            readerIndex++;
            return hasNext();

        }
        return nextPosRecord != null;
    }



    @Override
    public void remove() {
        throw new UnsupportedOperationException("Remove is not supported by this iterator.");
    }
}
