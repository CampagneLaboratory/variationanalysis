package org.campagnelab.dl.varanalysis.learning.iterators;

import it.unimi.dsi.util.XoRoShiRo128PlusRandom;
import org.nd4j.linalg.api.ndarray.INDArray;
import org.nd4j.linalg.dataset.DataSet;
import org.nd4j.linalg.dataset.api.DataSetPreProcessor;
import org.nd4j.linalg.dataset.api.iterator.DataSetIterator;
import org.nd4j.linalg.factory.Nd4j;

import java.io.Serializable;
import java.util.Arrays;
import java.util.Iterator;
import java.util.List;
import java.util.function.Consumer;

/**
 * An iterator that samples input records according to some probability.
 * Created by fac2003 on 7/21/16.
 */
public class SamplingIterator implements Iterator<DataSet>, org.nd4j.linalg.dataset.api.iterator.DataSetIterator, Serializable {
    private final float[] samplingProbabilities;
    private int miniBatchSize;
    private float averageProbability;
    private org.nd4j.linalg.dataset.api.iterator.DataSetIterator delegate;
    /**
     * Number of examples in the delegate iterator (before sampling).
     */
    private int numExamples;
    /**
     * The offset in samplingProbabilities of the first record of the last minibatch returned by next(int)
     */
    private int offsetStartOfMinibatch;
    /**
     * The size of the last minibatch returned by next().
     */
    private int lastMinibatchSize;
    private int delegateBatchSize;

    public SamplingIterator(DataSetIterator delegate, long seed) {
        this.delegate = delegate;
        this.randomGenerator = new XoRoShiRo128PlusRandom(seed);
        this.samplingProbabilities = new float[delegate.numExamples()];
        // sampling probabilities are set to 1, forcing all records to be returned by default.
        Arrays.fill(samplingProbabilities, 1f);
        recalculateMinibatchSize();
    }

    private void recalculateMinibatchSize() {

        miniBatchSize=delegate.batch();
        averageProbability = 0f;
        for (float p : samplingProbabilities) {
            averageProbability += p;
        }
        averageProbability /= samplingProbabilities.length;
        float willIgnoreFraction = 1 - averageProbability;
        this.delegateBatchSize = Math.round(delegate.batch() *(1+willIgnoreFraction));
        assert delegateBatchSize >= miniBatchSize : "delegate minibatch size must be larger than output size.";

    }

    public void setSamplingProbability(int indexInMinibatch, float probability) {
        samplingProbabilities[offsetStartOfMinibatch - lastMinibatchSize + indexInMinibatch] = probability;
    }

    @Override
    public boolean hasNext() {
        return delegate.hasNext();
    }

    @Override
    public DataSet next() {

        DataSet dataset = sample(delegate.next(delegateBatchSize), samplingProbabilities, offsetStartOfMinibatch);
        this.offsetStartOfMinibatch += dataset.numExamples();
        this.lastMinibatchSize = dataset.numExamples();
        return dataset;
    }

    private XoRoShiRo128PlusRandom randomGenerator;

    /**
     * Sample from delegate according to the sampling probabilities. Include a record from the delegate
     * only when a random number is less than the sampling probability.
     *
     * @param next
     * @param samplingProbabilities
     * @param offsetStartOfMinibatch
     * @return
     */
    private DataSet sample(DataSet next, float[] samplingProbabilities, int offsetStartOfMinibatch) {
        int numSamples = miniBatchSize;
        INDArray examples = Nd4j.create(numSamples, next.getFeatures().columns());
        INDArray outcomes = Nd4j.create(numSamples, next.numOutcomes());
        int delegateNumExamples = next.numExamples();
        int num = 0;
        int delegateIndex = 0;
        do {
            if (randomGenerator.nextFloat() < samplingProbabilities[offsetStartOfMinibatch + delegateIndex]) {

                // include in sample:
                DataSet dataSet = next.get(delegateIndex);
                examples.putRow(num, dataSet.getFeatures());
                outcomes.putRow(num, dataSet.getLabels());
                num++;
            }
            delegateIndex++;
            if (delegateIndex >= delegateNumExamples) break;
        } while (num < miniBatchSize );
        return new DataSet(examples, outcomes);
    }

    @Override
    public void remove() {
        delegate.remove();
    }

    @Override
    public void forEachRemaining(Consumer<? super DataSet> action) {
        delegate.forEachRemaining(action);
    }

    @Override
    public DataSet next(int num) {
        return delegate.next(miniBatchSize);
    }

    @Override
    public int totalExamples() {
        return delegate.totalExamples();
    }

    @Override
    public int inputColumns() {
        return delegate.inputColumns();
    }

    @Override
    public int totalOutcomes() {
        return delegate.totalOutcomes();
    }

    @Override
    public void reset() {
        recalculateMinibatchSize();
        delegate.reset();
    }

    @Override
    public int batch() {
        return delegate.batch();
    }

    @Override
    public int cursor() {
        return delegate.cursor();
    }

    @Override
    public int numExamples() {
        this.numExamples = delegate.numExamples();
        return numExamples;
    }

    @Override
    public void setPreProcessor(DataSetPreProcessor preProcessor) {
        delegate.setPreProcessor(preProcessor);
    }

    @Override
    public DataSetPreProcessor getPreProcessor() {
        return delegate.getPreProcessor();
    }

    @Override
    public List<String> getLabels() {
        return delegate.getLabels();
    }
}
