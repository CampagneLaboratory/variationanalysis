package org.campagnelab.dl.varanalysis.learning.iterators;

import org.apache.commons.compress.utils.IOUtils;
import org.campagnelab.dl.model.utils.mappers.FeatureMapper;
import org.campagnelab.dl.model.utils.mappers.LabelMapper;
import org.campagnelab.dl.varanalysis.protobuf.BaseInformationRecords;
import org.campagnelab.dl.varanalysis.storage.RecordReader;
import org.deeplearning4j.datasets.iterator.DataSetIterator;
import org.nd4j.linalg.api.ndarray.INDArray;
import org.nd4j.linalg.dataset.DataSet;
import org.nd4j.linalg.dataset.api.DataSetPreProcessor;
import org.nd4j.linalg.factory.Nd4j;

import java.io.IOException;
import java.util.List;
import java.util.NoSuchElementException;

/**
 * An iterator compatible with the DeepLearning4J framework that can iterate a BaseInformation file
 * and generate features suitable to train the neural net. The iterator can be constructed with a
 * dataset filename and an instance of FeatureCalculator.
 * <p>
 * Created by fac2003 on 5/21/16.
 *
 * @author Fabien Campagne
 */
public class BaseInformationIterator implements DataSetIterator {
    private final long totalExamples;
    private RecordReader reader;
    private final FeatureMapper featureMapper;
    private final LabelMapper labelMapper;
    private final String inputFilename;
    private BaseInformationRecords.BaseInformationOrBuilder nextPosRecord;
    private long cursor = 0;
    private int batchSize = 32;

    public BaseInformationIterator(String inputFilename, int batchSize, FeatureMapper featureMapper, LabelMapper labelMapper) throws IOException {
        this.inputFilename = inputFilename;
        this.featureMapper = featureMapper;
        this.labelMapper = labelMapper;
        this.batchSize = batchSize;
        this.reader = new RecordReader(inputFilename);
        this.totalExamples = reader.getTotalRecords();

    }

    @Override
    public DataSet next(int batchSize) {

        // allocate a new dataset with batchSize records and fill it with features and labels.

        int size = Math.min(batchSize, (int) remainingExamples());

        // allocate features and labels for the entire dataset:
        // dimension 0 = number of examples in minibatch
        // dimension 1 = number of features per record.

        //size changed from batchSize. huge batchSize values useful for tests
        INDArray inputs = Nd4j.zeros(size, featureMapper.numberOfFeatures());
        INDArray labels = Nd4j.zeros(size, labelMapper.numberOfLabels());
        for (int i = 0; i < size; i++) {
            // we are going to call nextRecord directly, without checking hasNextRecord, because we have
            // determined how many times we can call (in size). We should get the exception if we were
            // wrong in our estimate of size.

            // fill in features and labels for a given record i:
            BaseInformationRecords.BaseInformationOrBuilder record = nextRecord();
            featureMapper.mapFeatures(record, inputs, i);
            labelMapper.mapLabels(record, labels, i);

        }
        return new DataSet(inputs, labels);
    }


    private long remainingExamples() {
        return totalExamples - cursor;
    }

    @Override
    public int totalExamples() {
        return (int) totalExamples;
    }

    @Override
    public int inputColumns() {
        return featureMapper.numberOfFeatures();
    }

    @Override
    public int totalOutcomes() {
        return labelMapper.numberOfLabels();
    }

    @Override
    public void reset() {
        if (this.reader != null) {
            IOUtils.closeQuietly(reader);
        }
        try {
            this.reader = new RecordReader(inputFilename);
        }catch (IOException e) {

            throw new RuntimeException(e);
        }
        cursor = 0;
        nextPosRecord = null;
        //    System.out.println("reset called");
    }

    @Override
    public int batch() {
        return batchSize;
    }

    @Override
    public int cursor() {

        return (int) cursor;
    }

    @Override
    public int numExamples() {
        return (int) totalExamples;
    }

    @Override
    public void setPreProcessor(DataSetPreProcessor dataSetPreProcessor) {

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
        try {
            this.nextPosRecord = reader.nextRecord();
        } catch (IOException e) {
            throw new RuntimeException(e);
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
        try {
            this.nextPosRecord = reader.nextRecord();
        } catch (IOException e) {
            throw new RuntimeException(e);
        }
        return nextPosRecord != null;
    }

    /**
     * Return the next single PosRecord.
     *
     * @return the next available record, or throws NoSuchElementException if there are no more records.
     */
    BaseInformationRecords.BaseInformationOrBuilder nextRecord() {
        if (hasNextRecord()) {
            BaseInformationRecords.BaseInformationOrBuilder tmp = nextPosRecord;
            // setting nextPosRecord will make hasNextRecord load the next record from the underlying reader.
            nextPosRecord = null;
            cursor += 1;
            return tmp;
        } else throw new NoSuchElementException();
    }

    @Override
    public void remove() {
        throw new UnsupportedOperationException("Remove is not supported by this iterator.");
    }


}
