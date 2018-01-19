package org.campagnelab.dl.framework.iterators;

import it.unimi.dsi.fastutil.objects.ObjectArrayList;
import it.unimi.dsi.fastutil.objects.ObjectList;
import org.campagnelab.dl.framework.domains.DomainDescriptor;
import org.campagnelab.dl.framework.mappers.ConfigurableFeatureMapper;
import org.campagnelab.dl.framework.mappers.FeatureMapper;
import org.campagnelab.dl.framework.mappers.LabelMapper;
import org.campagnelab.goby.baseinfo.SequenceBaseInformationReader;
import org.nd4j.linalg.api.ndarray.INDArray;
import org.nd4j.linalg.dataset.api.MultiDataSet;
import org.nd4j.linalg.dataset.api.MultiDataSetPreProcessor;
import org.nd4j.linalg.dataset.api.iterator.MultiDataSetIterator;
import org.nd4j.linalg.factory.Nd4j;

import java.io.IOException;
import java.lang.reflect.Constructor;
import java.lang.reflect.InvocationTargetException;
import java.lang.reflect.Type;
import java.util.*;

/**
 * Make a multi dataset iterator from an iterable over records.
 */
public abstract class MultiDataSetIteratorAdapterMultipleSamples<RecordType> implements Iterator<List<MultiDataSet>>,
        Iterable<List<MultiDataSet>> {

    private final DomainDescriptor domainDescriptor;
    private final Iterable<RecordType> iterable;
    private final int[] sampleIndices;
    private Iterator<RecordType> recordIterator;
    private boolean isPretrained;
    private Integer eosIndex;
    private FeatureMapper[][] featureMappers;
    private LabelMapper[][] labelMappers;

    protected long totalExamples;


    protected int batchSize = 32;
    private MultiDataSetPreProcessor preProcessor;


    public MultiDataSetIteratorAdapterMultipleSamples(Iterable<RecordType> iterable, int batchSize, DomainDescriptor domainDescriptor,
                                                      boolean isPretrained, Integer eosIndex, int[] sampleIndices,
                                                      Properties readerProperties) throws IOException {
        this.domainDescriptor = domainDescriptor;
        this.batchSize = batchSize;
        this.iterable = iterable;
        this.recordIterator = iterable.iterator();
        this.isPretrained = isPretrained;
        this.eosIndex = eosIndex;
        this.sampleIndices = sampleIndices;
        final int numInputs = domainDescriptor.getComputationalGraph().getInputNames().length;
        final int numLabels = domainDescriptor.getComputationalGraph().getOutputNames().length;
        featureMappers = new FeatureMapper[numInputs][sampleIndices.length];
        labelMappers = new LabelMapper[numLabels][sampleIndices.length];
        int index = 0;
        for (String input : domainDescriptor.getComputationalGraph().getInputNames()) {
            int sampleIdIndex = 0;
            for (int sampleIndex : sampleIndices) {

                FeatureMapper sampleFeatureMapper = domainDescriptor.getFeatureMapper(input, sampleIndex);
                if (sampleFeatureMapper instanceof ConfigurableFeatureMapper) {
                    ConfigurableFeatureMapper cmapper = (ConfigurableFeatureMapper) sampleFeatureMapper;
                    cmapper.configure(readerProperties);
                }
                featureMappers[index][sampleIdIndex] = sampleFeatureMapper;

            }
            index++;
        }
        index = 0;
        for (String label : domainDescriptor.getComputationalGraph().getOutputNames()) {
            for (int sampleIndex : sampleIndices) {
                labelMappers[index][sampleIndex] = domainDescriptor.getLabelMapper(label, sampleIndex);
            }
            index++;
        }
    }

    abstract public String getBasename();


    public List<MultiDataSet> next(int batchSize) {
        if (!hasNext()) {
            throw new NoSuchElementException();
        }
        ObjectList<RecordType> buffer = new ObjectArrayList<RecordType>();
        // allocate a new dataset with batchSize records and fill it with features and labels.
        while (recordIterator.hasNext() && buffer.size() < batchSize) {
            buffer.add(recordIterator.next());
        }
        int size = buffer.size();

        // allocate features and labels for the entire dataset:
        // dimension 0 = number of examples in minibatch
        // dimension 1 = number of features per record.
        final int numInputs = domainDescriptor.getComputationalGraph().getInputNames().length;
        final int numLabels = domainDescriptor.getComputationalGraph().getOutputNames().length;

        INDArray inputs[][] = new INDArray[sampleIndices.length][numInputs];
        INDArray inputMasks[][] = new INDArray[sampleIndices.length][numInputs];
        INDArray labels[][] = new INDArray[sampleIndices.length][numLabels];
        INDArray labelMasks[][] = new INDArray[sampleIndices.length][numLabels];
        int index = 0;
        boolean hasFeatureMask = false;
        boolean hasLabelMask = false;

        for (String input : domainDescriptor.getComputationalGraph().getInputNames()) {
            int[] inputShape = domainDescriptor.getInputShape(size, input).clone();
            boolean padEos = (isPretrained) && ((eosIndex != null && eosIndex == inputShape[1]) || eosIndex == null);
            if (padEos) {
                if (inputShape.length != 3) {
                    throw new RuntimeException("EOS padding only valid for sequences with 2D features");
                }
                inputShape[1]++;
            }
            boolean needMask = featureMappers[index][0].hasMask();
            for (int sampleIndex = 0; sampleIndex < sampleIndices.length; sampleIndex++) {
                inputs[sampleIndex][index] = Nd4j.createUninitializedDetached(inputShape, 'f');
                inputMasks[sampleIndex][index] = needMask ? Nd4j.createUninitializedDetached(
                        domainDescriptor.getInputMaskShape(size, input), 'f'
                ) : null;
            }
            index++;
            hasFeatureMask |= needMask;
        }
        index = 0;
        for (String label : domainDescriptor.getComputationalGraph().getOutputNames()) {
            boolean needMask = labelMappers[index][0].hasMask();
            for (int sampleIndex = 0; sampleIndex < sampleIndices.length; sampleIndex++) {
                labelMappers[index][sampleIndex] = domainDescriptor.getLabelMapper(label);
                labels[sampleIndex][index] = Nd4j.createUninitializedDetached(domainDescriptor.getLabelShape(size, label), 'f');
                labelMasks[sampleIndex][index] = needMask ? Nd4j.createUninitializedDetached(
                        domainDescriptor.getLabelMaskShape(size, label), 'f'
                ) : null;
            }
            index++;
            hasLabelMask |= needMask;
        }


        int recordIndexInBatch = 0;
        for (RecordType record : buffer) {
            for (int j = 0; j < numInputs; j++) {
                for (int sampleIndex = 0; sampleIndex < sampleIndices.length; sampleIndex++) {
                    featureMappers[j][sampleIndex].prepareToNormalize(record, recordIndexInBatch);
                    featureMappers[j][sampleIndex].mapFeatures(record, inputs[sampleIndex][j], recordIndexInBatch);
                    if (featureMappers[j][sampleIndex].hasMask()) {
                        featureMappers[j][sampleIndex].maskFeatures(record, inputMasks[sampleIndex][j],
                                recordIndexInBatch);
                    }
                }
            }
            for (int j = 0; j < numLabels; j++) {
                for (int sampleIndex = 0; sampleIndex < sampleIndices.length; sampleIndex++) {
                    labelMappers[j][sampleIndex].prepareToNormalize(record, recordIndexInBatch);
                    labelMappers[j][sampleIndex].mapLabels(record, labels[sampleIndex][j], recordIndexInBatch);
                    if (labelMappers[j][sampleIndex].hasMask()) {
                        labelMappers[j][sampleIndex].maskLabels(record, labelMasks[sampleIndex][j],
                                recordIndexInBatch);
                    }
                }
            }
            recordIndexInBatch += 1;
        }
        // Necessary for mixed datasets (i.e., where some mappers have masks and others don't) - will raise NPE otherwise
        if (hasFeatureMask) {
            for (int i = 0; i < inputMasks[0].length; i++) {
                if (inputMasks[0][i] == null) {
                    int[] inputShape = inputs[0][i].shape();
                    if (inputShape.length == 3) {
                        //     inputMasks[i] = Nd4j.ones(inputShape[0], 1,inputShape[2]);
                        throw new RuntimeException("3D features should have masks");
                    } else if (inputShape.length == 2 || inputShape.length == 1) {
                        for (int sampleIndex = 0; sampleIndex < sampleIndices.length; sampleIndex++)
                            inputMasks[sampleIndex][i] = Nd4j.ones(inputShape[0], 1);
                    } else {
                        for (int sampleIndex = 0; sampleIndex < sampleIndices.length; sampleIndex++)
                            inputMasks[sampleIndex][i] = Nd4j.ones(inputShape.clone());
                    }
                }
            }
        }
        if (hasLabelMask) {
            for (int i = 0; i < labelMasks.length; i++) {
                if (labelMasks[i] == null) {
                    int[] labelShape = labels[0][i].shape();
                    if (labelShape.length == 3) {
                        //  labelMasks[i] = Nd4j.ones(labelShape[0], labelShape[2], 1);
                        //throw new RuntimeException("3D labels should have masks");
                    } else if (labelShape.length == 2 || labelShape.length == 1) {
                        for (int sampleIndex = 0; sampleIndex < sampleIndices.length; sampleIndex++)
                            labelMasks[sampleIndex][i] = Nd4j.ones(labelShape[0], 1);
                    } else {
                        for (int sampleIndex = 0; sampleIndex < sampleIndices.length; sampleIndex++)
                            labelMasks[sampleIndex][i] = Nd4j.ones(labelShape.clone());
                    }
                }
            }
        }
        List<MultiDataSet> resultList = new ArrayList<>();
        for (int sampleIndex = 0; sampleIndex < sampleIndices.length; sampleIndex++) {
            MultiDataSet result = new org.nd4j.linalg.dataset.MultiDataSet(inputs[sampleIndex], labels[sampleIndex],
                    hasFeatureMask ? inputMasks[sampleIndex] : null,
                    hasLabelMask ? labelMasks[sampleIndex] : null);
            if (preProcessor != null) preProcessor.preProcess(result);
            resultList.add(result);
        }
        return resultList;
    }


    public boolean asyncSupported() {
        return true;
    }

    public void reset() {
        recordIterator = iterable.iterator();
    }

    public void setPreProcessor(MultiDataSetPreProcessor preProcessor) {
        this.preProcessor = preProcessor;
    }

    public MultiDataSetPreProcessor getPreProcessor() {
        return this.preProcessor;
    }

    @Override
    public boolean hasNext() {
        return recordIterator.hasNext();
    }

    @Override
    public List<MultiDataSet> next() {
        if (hasNext()) {
            return next(batchSize);

        } else throw new NoSuchElementException();
    }

    @Override
    public void remove() {
        throw new UnsupportedOperationException("Remove is not supported by this iterator.");
    }

    @Override
    public Iterator<List<MultiDataSet>> iterator() {
        reset();
        return this;
    }
}