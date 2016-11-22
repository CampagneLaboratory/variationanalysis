package org.campagnelab.dl.framework.iterators;

import it.unimi.dsi.fastutil.objects.ObjectArrayList;
import it.unimi.dsi.fastutil.objects.ObjectList;
import org.campagnelab.dl.framework.domains.DomainDescriptor;
import org.campagnelab.dl.framework.mappers.FeatureMapper;
import org.campagnelab.dl.framework.mappers.LabelMapper;
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
public abstract class MultiDataSetIteratorAdapter<RecordType> implements MultiDataSetIterator, Iterable<MultiDataSet> {

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

        INDArray inputs[] = new INDArray[numInputs];
        INDArray inputMasks[] = new INDArray[numInputs];
        INDArray labels[] = new INDArray[numLabels];
        INDArray labelMasks[] = new INDArray[numLabels];
        FeatureMapper[] featureMappers = new FeatureMapper[numInputs];
        LabelMapper[] labelMappers = new LabelMapper[numOutputs];
        int index = 0;
        boolean hasFeatureMask = false;
        boolean hasLabelMask = false;
        for (String input : domainDescriptor.getComputationalGraph().getInputNames()) {

            inputs[index] = Nd4j.zeros(domainDescriptor.getInputShape(size, input));
            featureMappers[index] = domainDescriptor.getFeatureMapper(input);
            boolean needMask = featureMappers[index].hasMask();
            inputMasks[index] = needMask ? Nd4j.zeros(domainDescriptor.getInputMaskShape(size, input)) : null;
            index += 1;
            hasFeatureMask |= needMask;
        }
        index = 0;
        for (String label : domainDescriptor.getComputationalGraph().getOutputNames()) {
            labels[index] = Nd4j.zeros(domainDescriptor.getLabelShape(size, label));
            labelMappers[index] = domainDescriptor.getLabelMapper(label);
            boolean needMask = labelMappers[index].hasMask();
            labelMasks[index] = needMask ? Nd4j.zeros(domainDescriptor.getLabelMaskShape(size, label)) : null;
            index++;
            hasLabelMask |= needMask;
        }
        int recordIndexInBatch = 0;
        for (RecordType record : buffer) {

            for (int j = 0; j < numInputs; j++) {
                featureMappers[j].prepareToNormalize(record, recordIndexInBatch);
                featureMappers[j].mapFeatures(record, inputs[j], recordIndexInBatch);
                if (featureMappers[j].hasMask()) {
                    featureMappers[j].maskFeatures(record, inputMasks[j], recordIndexInBatch);
                }
            }
            for (int j = 0; j < numOutputs; j++) {
                labelMappers[j].mapLabels(record, labels[j], recordIndexInBatch);
                if (labelMappers[j].hasMask()) {
                    labelMappers[j].maskLabels(record, labelMasks[j], recordIndexInBatch);
                }
            }
            recordIndexInBatch += 1;

        }
        final MultiDataSet result = new org.nd4j.linalg.dataset.MultiDataSet(inputs, labels,
                hasFeatureMask ? inputMasks : null,
                hasLabelMask ? labelMasks : null);
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

    @Override
    public Iterator<MultiDataSet> iterator() {
        reset();
        return this;
    }
}
