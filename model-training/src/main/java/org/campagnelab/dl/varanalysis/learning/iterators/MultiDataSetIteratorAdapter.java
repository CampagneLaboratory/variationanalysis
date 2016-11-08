package org.campagnelab.dl.varanalysis.learning.iterators;

import it.unimi.dsi.fastutil.objects.ObjectArrayList;
import it.unimi.dsi.fastutil.objects.ObjectList;
import org.campagnelab.dl.model.utils.mappers.FeatureMapper;
import org.campagnelab.dl.model.utils.mappers.LabelMapper;
import org.campagnelab.dl.varanalysis.learning.DomainDescriptor;
import org.nd4j.linalg.api.ndarray.INDArray;
import org.nd4j.linalg.dataset.api.MultiDataSet;
import org.nd4j.linalg.dataset.api.MultiDataSetPreProcessor;
import org.nd4j.linalg.dataset.api.iterator.MultiDataSetIterator;
import org.nd4j.linalg.factory.Nd4j;

import java.io.IOException;
import java.util.Iterator;
import java.util.NoSuchElementException;

/**
 * Make a multi dataset iterator from an iterable over records.
 */
public abstract class MultiDataSetIteratorAdapter<RecordType> implements MultiDataSetIterator {

    private final DomainDescriptor domainDescriptor;
    private final Iterable<RecordType> iterable;
    private Iterator<RecordType> recordIterator;

    protected long totalExamples;


    protected int batchSize = 32;
    private MultiDataSetPreProcessor preProcessor;


    public MultiDataSetIteratorAdapter(Iterable<RecordType> iterable, int batchSize, DomainDescriptor domainDescriptor) throws IOException {
        this.domainDescriptor = domainDescriptor;
        this.batchSize = batchSize;
        this.iterable = iterable;
        this.recordIterator = iterable.iterator();
    }

    abstract public String getBasename();


    ObjectList<RecordType> buffer = new ObjectArrayList<RecordType>();

    public MultiDataSet next(int batchSize) {
        buffer.clear();
        // allocate a new dataset with batchSize records and fill it with features and labels.
        while (recordIterator.hasNext() && buffer.size() < batchSize) {
            buffer.add(recordIterator.next());
        }
        int size = buffer.size();

        // allocate features and labels for the entire dataset:
        // dimension 0 = number of examples in minibatch
        // dimension 1 = number of features per record.

        //size changed from batchSize. huge batchSize values useful for tests
        final int numInputs = domainDescriptor.getComputationalGraph().getInputNames().length;
        final int numLabels = domainDescriptor.getComputationalGraph().getOutputNames().length;
        int numOutputs = numLabels;

        INDArray inputs[] = new INDArray[numInputs];//= Nd4j.zeros(size, featureMapper.numberOfFeatures());
        INDArray labels[] = new INDArray[numLabels];//= Nd4j.zeros(size, labelMapper.numberOfLabels());
        FeatureMapper[] featureMappers = new FeatureMapper[numInputs];
        LabelMapper[] labelMappers = new LabelMapper[numOutputs];
        int index = 0;

        for (String input : domainDescriptor.getComputationalGraph().getInputNames()) {

            inputs[index] = Nd4j.zeros(domainDescriptor.getInputShape(size, input));
            featureMappers[index] = domainDescriptor.getFeatureMapper(input);
            index += 1;

        }
        index = 0;
        for (String label : domainDescriptor.getComputationalGraph().getOutputNames()) {
            labels[index] = Nd4j.zeros(domainDescriptor.getLabelShape(size, label));
            labelMappers[index] = domainDescriptor.getLabelMapper(label);
            index++;
        }
        int recordIndexInBatch = 0;
        for (RecordType record : buffer) {

            for (int j = 0; j < numInputs; j++) {
                featureMappers[j].prepareToNormalize(record, recordIndexInBatch);
                featureMappers[j].mapFeatures(record, inputs[j], recordIndexInBatch);
            }
            for (int j = 0; j < numOutputs; j++) {
                labelMappers[j].mapLabels(record, labels[j], recordIndexInBatch);
            }
            recordIndexInBatch += 1;

        }
        final MultiDataSet result = new org.nd4j.linalg.dataset.MultiDataSet(inputs, labels);
        if (preProcessor != null) preProcessor.preProcess(result);
        return result;
    }

    @Override
    public void setPreProcessor(MultiDataSetPreProcessor preProcessor) {
        this.preProcessor = preProcessor;
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

        recordIterator = iterable.iterator();
    }


    @Override
    public boolean hasNext() {
        return recordIterator.hasNext();
    }


    @Override
    public MultiDataSet next() {
        if (hasNext()) {
            return next(batchSize);

        } else throw new NoSuchElementException();
    }

    @Override
    public void remove() {
        throw new UnsupportedOperationException("Remove is not supported by this iterator.");
    }

}
