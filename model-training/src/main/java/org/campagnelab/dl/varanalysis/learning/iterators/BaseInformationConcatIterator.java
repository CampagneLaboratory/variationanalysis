package org.campagnelab.dl.varanalysis.learning.iterators;


import it.unimi.dsi.fastutil.objects.ObjectArrayList;
import org.campagnelab.dl.model.utils.mappers.FeatureMapper;
import org.campagnelab.dl.model.utils.mappers.LabelMapper;
import org.nd4j.linalg.dataset.DataSet;
import org.nd4j.linalg.dataset.api.DataSetPreProcessor;

import java.io.IOException;
import java.util.List;
import java.util.NoSuchElementException;

/**
 * Created by rct66 on 6/13/16.
 */
public class BaseInformationConcatIterator extends BaseInformationIterator implements NamedDataSetIterator {

    private List<BaseInformationIterator> baseIters;
    private int readerIndex = 0;

    public BaseInformationConcatIterator(List<BaseInformationIterator> iterators, int batchSize, FeatureMapper featureMapper, LabelMapper labelMapper) throws IOException {
        super(featureMapper, labelMapper);
        this.batchSize = batchSize;
        this.baseIters = iterators;
        for (BaseInformationIterator iter : iterators) {
            this.totalExamples += iter.totalExamples();
        }

    }

    public BaseInformationConcatIterator(int batchSize, FeatureMapper featureMapper, LabelMapper labelMapper, String... inputFilename) throws IOException {
        this(buildIterators(inputFilename, batchSize, featureMapper, labelMapper), batchSize, featureMapper, labelMapper);
    }

    private static List<BaseInformationIterator> buildIterators(String[] inputFilename, int batchSize, FeatureMapper featureMapper, LabelMapper labelMapper) throws IOException {
        final ObjectArrayList list = new ObjectArrayList();
        for (String filename : inputFilename) {

            list.add(new BaseInformationIterator(filename, batchSize, featureMapper, labelMapper));
        }
        return list;
    }

    @Override
    public String getBasename() {
        if (baseIters.size() == 1) {
            return baseIters.get(0).getBasename();
        } else {
            String basename = null;
            long hashcode = 8723872838723L;
            for (NamedDataSetIterator it : baseIters) {
                hashcode ^= it.getBasename().hashCode();
            }
            basename = "multiset-" + Long.toString(hashcode);
            return basename;
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
        if (readerIndex >= baseIters.size()) {
            return false;
        }
        try {
            this.nextPosRecord = baseIters.get(readerIndex).nextRecord();
        } catch (NoSuchElementException e) {
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
        if (readerIndex >= baseIters.size()) {
            return false;
        }
        try {
            this.nextPosRecord = baseIters.get(readerIndex).nextRecord();
        } catch (NoSuchElementException e) {
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
