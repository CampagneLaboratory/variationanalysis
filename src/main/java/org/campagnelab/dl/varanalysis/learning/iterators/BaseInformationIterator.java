package org.campagnelab.dl.varanalysis.learning.iterators;

import org.apache.commons.compress.utils.IOUtils;
import org.campagnelab.dl.varanalysis.format.PosRecord;
import org.campagnelab.dl.varanalysis.storage.AvroVariationParquetReader;
import org.deeplearning4j.datasets.iterator.DataSetIterator;
import org.nd4j.linalg.api.ndarray.INDArray;
import org.nd4j.linalg.dataset.DataSet;
import org.nd4j.linalg.dataset.api.DataSetPreProcessor;
import org.nd4j.linalg.factory.Nd4j;

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
    private final int totalExamples;
    private AvroVariationParquetReader reader;
    private final FeatureCalculator featureCalculator;
    private final String inputFilename;
    private PosRecord nextPosRecord;
    private int cursor = 0;
    private int batchSize = 32;

    public BaseInformationIterator(String inputFilename, int batchSize, FeatureCalculator featureCalculator) {
        this.inputFilename = inputFilename;
        this.featureCalculator = featureCalculator;
        this.batchSize = batchSize;
        AvroVariationParquetReader reader = new AvroVariationParquetReader(inputFilename);
        // traverse the data once to find the total number of records:
        int count = 0;
        try {

            while
                    (reader.read() != null) {
                count++;
            }

        } finally {
            reader.close();
        }
        this.totalExamples = count;
        this.reader = new AvroVariationParquetReader(inputFilename);
    }

    @Override
    public DataSet next(int batchSize) {

        // allocate a new dataset with batchSize records and fill it with features and labels.
        DataSet ds = new DataSet();
        int size = Math.min(batchSize, remainingExamples());

        // allocate features and labels for the entire dataset:
        // dimension 0 = number of examples in minibatch
        // dimension 1 = number of features per record.

        INDArray inputs = Nd4j.zeros(batchSize, featureCalculator.numberOfFeatures());
        INDArray labels = Nd4j.zeros(batchSize, featureCalculator.numberOfLabels());
        for (int i = 0; i < size; i++) {
            // we are going to call nextRecord directly, without checking hasNextRecord, because we have
            // determined how many times we can call (in size). We should get the exception if we were
            // wrong in our estimate of size.

            // fill in features and labels for a given record i:
            PosRecord record = nextRecord();
            featureCalculator.map(record, inputs, labels, i);
            ds.setFeatures(inputs);
            ds.setLabels(labels);
        }
        return ds;
    }


    private int remainingExamples() {
        return totalExamples - cursor;
    }

    @Override
    public int totalExamples() {
        return totalExamples;
    }

    @Override
    public int inputColumns() {
        return featureCalculator.numberOfFeatures();
    }

    @Override
    public int totalOutcomes() {
        return featureCalculator.numberOfLabels();
    }

    @Override
    public void reset() {
        if (this.reader != null) {
            IOUtils.closeQuietly(reader);
        }
        this.reader = new AvroVariationParquetReader(inputFilename);
        cursor = 0;
    }

    @Override
    public int batch() {
        return batchSize;
    }

    @Override
    public int cursor() {

        return cursor;
    }

    @Override
    public int numExamples() {
        return totalExamples;
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
        this.nextPosRecord = reader.read();
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
        this.nextPosRecord = reader.read();
        return nextPosRecord != null;
    }

    /**
     * Return the next single PosRecord.
     *
     * @return the next available record, or throws NoSuchElementException if there are no more records.
     */
    private PosRecord nextRecord() {
        if (hasNext()) {
            PosRecord tmp = nextPosRecord;
            // setting nextPosRecord will make hasNextRecord load the next record from the underlying reader.
            nextPosRecord = null;
            return tmp;
        } else throw new NoSuchElementException();
    }

    @Override
    public void remove() {
        throw new UnsupportedOperationException("Remove is not supported by this iterator.");
    }


}
