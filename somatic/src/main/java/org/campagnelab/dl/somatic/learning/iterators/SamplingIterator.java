package org.campagnelab.dl.somatic.learning.iterators;

import it.unimi.dsi.fastutil.ints.Int2IntArrayMap;
import it.unimi.dsi.fastutil.ints.IntArrayList;
import it.unimi.dsi.util.XoRoShiRo128PlusRandom;
import org.nd4j.linalg.api.ndarray.INDArray;
import org.nd4j.linalg.dataset.DataSet;
import org.nd4j.linalg.dataset.api.DataSetPreProcessor;
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
public class SamplingIterator implements Iterator<DataSet>, NamedDataSetIterator, Serializable {
    private final float[] samplingProbabilities;
    private int miniBatchSize;
    private float averageProbability;
    private NamedDataSetIterator delegate;
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
    private int numReturned;
    private int numSkipped;
    private int currentRecordIndex;

    public SamplingIterator(NamedDataSetIterator delegate, long seed) {

        this.delegate = delegate;
        this.randomGenerator = new XoRoShiRo128PlusRandom(seed);
        this.samplingProbabilities = new float[delegate.numExamples()];
        // sampling probabilities are set to 1, forcing all records to be returned by default.
        Arrays.fill(samplingProbabilities, 1f);
        recalculateMinibatchSize();
        throw new UnsupportedOperationException("This iterator is broken in DL4J 0.6.0, do not use.");
    }

    public void updateStatistics() {
        recalculateMinibatchSize();
    }

    private void recalculateMinibatchSize() {

        miniBatchSize = delegate.batch();
        averageProbability = 0f;
        for (float p : samplingProbabilities) {
            averageProbability += p;
        }
        averageProbability /= samplingProbabilities.length;
        float willIgnoreFraction = 1 - averageProbability;
        this.delegateBatchSize = Math.round(delegate.batch() * (1 / averageProbability));
        assert delegateBatchSize >= miniBatchSize : "delegate minibatch size must be larger than output size.";

    }

    public float getAverageSamplingP() {
        return averageProbability;
    }

    public void setSamplingProbability(boolean wrongPrediction, int indexInMinibatch, float probability) {
        if (!sampleToDelegateIndexMap.containsKey(indexInMinibatch)) {
            System.out.println("STOP");
        }
        assert sampleToDelegateIndexMap.containsKey(indexInMinibatch) : "The index (in minibatch) must be found: " + indexInMinibatch;
        final int index = sampleToDelegateIndexMap.get(indexInMinibatch);
        if (index >= samplingProbabilities.length) {

        } else {
//            if (wrongPrediction) {
//                System.out.printf("Wrong prediction getting sample prob=%f %n",probability);
//            }
            samplingProbabilities[index] = probability;
        }

    }

    @Override
    public boolean hasNext() {
        return delegate.hasNext();
    }

    @Override
    public DataSet next() {
        //int previousPoint = numReturned + numSkipped;
        DataSet dataset = sample(delegate.next(delegateBatchSize), samplingProbabilities, offsetStartOfMinibatch);
        //int newPoint = numReturned + numSkipped;
        this.offsetStartOfMinibatch += lastMinibatchSize;
        //    this.lastMinibatchSize = delegateBatchSize;
        return dataset;
    }

    private XoRoShiRo128PlusRandom randomGenerator;
    private Int2IntArrayMap sampleToDelegateIndexMap = new Int2IntArrayMap();

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
        sampleToDelegateIndexMap.clear();
        int numSamples = miniBatchSize;
        lastMinibatchSize = delegateBatchSize;

        int delegateNumExamples = next.numExamples();
        lastMinibatchSize = delegateNumExamples;
        int num = 0;
        IntArrayList selectedIndices = new IntArrayList();
        for (int delegateIndex = 0; delegateIndex < delegateNumExamples; delegateIndex++) {

            final float random = randomGenerator.nextFloat();
            final float p = getSamplingProbability(currentRecordIndex);
            //   System.out.printf("random=%f p=%f%n",random,p);
            //   System.out.flush();

            if (random < p) {
                selectedIndices.add(delegateIndex);
                sampleToDelegateIndexMap.put(num, currentRecordIndex);
                numReturned++;
                num++;
            } else {
                numSkipped++;
            }
            currentRecordIndex++;

        }

        numSamples = selectedIndices.size();
        if (numSamples == 0) {
            return next;
        }
        num = 0;
        INDArray examples = Nd4j.create(numSamples, next.getFeatures().columns());
        INDArray outcomes = Nd4j.create(numSamples, next.numOutcomes());
        for (int delegateIndex : selectedIndices) {
            DataSet dataSet = next.get(delegateIndex);
            examples.putRow(num, dataSet.getFeatures());
            outcomes.putRow(num, dataSet.getLabels());
            num++;
        }

        return new DataSet(examples, outcomes);
    }

    public float getProbability(int indexInMinibatch) {
        assert sampleToDelegateIndexMap.containsKey(indexInMinibatch) : "The index (in mini-batch returned) must be found: " + indexInMinibatch;
        final int index = sampleToDelegateIndexMap.get(indexInMinibatch);
        if (index >= samplingProbabilities.length) {
            return 1;
        } else {
            return samplingProbabilities[index];
        }
    }

    private float getSamplingProbability(int i) {
        if (i >= samplingProbabilities.length) {
            //     System.out.println("STOP " + i);
            return 1;
        } else {
            return samplingProbabilities[i];
        }
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
    public boolean resetSupported() {
        return true;
    }


    public boolean asyncSupported() {
        return false;
    }

    @Override
    public void reset() {
        recalculateMinibatchSize();
        delegate.reset();
        numReturned = 0;
        numSkipped = 0;
        offsetStartOfMinibatch = 0;
        lastMinibatchSize = 0;
        currentRecordIndex = 0;
        sampleToDelegateIndexMap.clear();
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


    public double percentSkipped() {
        final double numSkipped = this.numSkipped;
        return 100f * numSkipped / (numSkipped + numReturned);
    }

    @Override
    public String getBasename() {
        return delegate.getBasename();
    }
}
