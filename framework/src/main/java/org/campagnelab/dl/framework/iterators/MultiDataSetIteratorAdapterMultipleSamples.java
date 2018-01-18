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
public abstract class MultiDataSetIteratorAdapterMultipleSamples<RecordType> implements Iterator<List<MultiDataSet>>, Iterable<List<MultiDataSet>> {

    private final DomainDescriptor domainDescriptor;
    private final Iterable<RecordType> iterable;
    private final Integer[] sampleIds;
    private Iterator<RecordType> recordIterator;
    private boolean isPretrained;
    private Integer eosIndex;
    private FeatureMapper[][] featureMappers;
    private LabelMapper[] labelMappers;

    protected long totalExamples;


    protected int batchSize = 32;
    private MultiDataSetPreProcessor preProcessor;


    public MultiDataSetIteratorAdapterMultipleSamples(Iterable<RecordType> iterable, int batchSize, DomainDescriptor domainDescriptor,
                                                      boolean isPretrained, Integer eosIndex, Integer[] sampleIds,
                                                      Properties readerProperties) throws IOException {
        this.domainDescriptor = domainDescriptor;
        this.batchSize = batchSize;
        this.iterable = iterable;
        this.recordIterator = iterable.iterator();
        this.isPretrained = isPretrained;
        this.eosIndex = eosIndex;
        this.sampleIds = sampleIds;
        final int numInputs = domainDescriptor.getComputationalGraph().getInputNames().length;
        final int numLabels = domainDescriptor.getComputationalGraph().getOutputNames().length;
        featureMappers = new FeatureMapper[numInputs][sampleIds.length];
        labelMappers = new LabelMapper[numLabels];
        int index = 0;
        for (String input : domainDescriptor.getComputationalGraph().getInputNames()) {
            FeatureMapper featureMapper = domainDescriptor.getFeatureMapper(input);
            Class featureMapperClass = featureMapper.getClass();
            Constructor[] featureMapperConstructors = featureMapperClass.getDeclaredConstructors();
            Constructor featureMapperSampleCtor = null;
            for (Constructor constructor : featureMapperConstructors) {
                if (constructor.getParameterCount() == 1) {
                    Class parameterType = constructor.getParameterTypes()[0];
                    if (parameterType == Integer.class) {
                        featureMapperSampleCtor = constructor;
                        break;
                    }
                }
            }
            if (featureMapperSampleCtor == null) {
                throw new IllegalArgumentException("Feature mapper class for input doesn't support sample id");
            }
            int sampleIdIndex = 0;
            for (int sampleId : sampleIds) {
                try {
                    FeatureMapper sampleFeatureMapper = (FeatureMapper) featureMapperSampleCtor.newInstance(sampleId);
                    if (featureMapper instanceof ConfigurableFeatureMapper) {
                        ConfigurableFeatureMapper cmapper = (ConfigurableFeatureMapper) sampleFeatureMapper;
                        cmapper.configure(readerProperties);
                    }
                    featureMappers[index][sampleIdIndex] = sampleFeatureMapper;
                } catch (InvocationTargetException e) {
                    throw new RuntimeException("Unable to instantiate or configure feature mapper", e);
                } catch (IllegalAccessException | InstantiationException e) {
                    throw new RuntimeException("IO excpetion, perhaps sbi file not found?", e);
                }
            }
            index++;
        }
        index = 0;
        for (String label : domainDescriptor.getComputationalGraph().getOutputNames()) {
            labelMappers[index] = domainDescriptor.getLabelMapper(label);
            index++;
        }
    }

    abstract public String getBasename();

    @Override
    public List<MultiDataSet> next() {
        if (!hasNext()) {
            throw new NoSuchElementException();
        }
        ObjectList<RecordType> buffer = new ObjectArrayList<RecordType>();
        // allocate a new dataset with batchSize records and fill it with features and labels.
        while (recordIterator.hasNext() && buffer.size() < this.batchSize) {
            buffer.add(recordIterator.next());
        }
        int size = buffer.size();

        // allocate features and labels for the entire dataset:
        // dimension 0 = number of examples in minibatch
        // dimension 1 = number of features per record.
        final int numInputs = domainDescriptor.getComputationalGraph().getInputNames().length;
        final int numLabels = domainDescriptor.getComputationalGraph().getOutputNames().length;

        INDArray inputs[][] = new INDArray[sampleIds.length][numInputs];
        INDArray inputMasks[][] = new INDArray[sampleIds.length][numInputs];
        INDArray labels[] = new INDArray[numLabels];
        INDArray labelMasks[] = new INDArray[numLabels];
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
            for (int i = 0; i < sampleIds.length; i++) {
                inputs[index][i] = Nd4j.createUninitializedDetached(inputShape, 'f');
                inputMasks[index][i] = needMask ? Nd4j.createUninitializedDetached(domainDescriptor.getInputMaskShape(size, input), 'f') : null;
            }


            index += 1;
            hasFeatureMask |= needMask;
        }
        index = 0;
        for (String label : domainDescriptor.getComputationalGraph().getOutputNames()) {
            labels[index] = Nd4j.createUninitializedDetached(domainDescriptor.getLabelShape(size, label), 'f');

            labelMappers[index] = domainDescriptor.getLabelMapper(label);
            boolean needMask = labelMappers[index].hasMask();
            if (needMask) {
                labelMasks[index] = Nd4j.createUninitializedDetached(domainDescriptor.getLabelMaskShape(size, label), 'f');
            }
            index++;
            hasLabelMask |= needMask;
        }

        int recordIndexInBatch = 0;
        for (RecordType record : buffer) {
            for (int j = 0; j < numInputs; j++) {
                for (int k = 0; k < sampleIds.length; k++) {
                    featureMappers[j][k].prepareToNormalize(record, recordIndexInBatch);
                    featureMappers[j][k].mapFeatures(record, inputs[k][j], recordIndexInBatch);
                    if (featureMappers[j][k].hasMask()) {
                        featureMappers[j][k].maskFeatures(record, inputMasks[k][j], recordIndexInBatch);
                    }
                }
            }
            for (int j = 0; j < numLabels; j++) {
                labelMappers[j].prepareToNormalize(record, recordIndexInBatch);
                labelMappers[j].mapLabels(record, labels[j], recordIndexInBatch);
                if (labelMappers[j].hasMask()) {
                    labelMappers[j].maskLabels(record, labelMasks[j], recordIndexInBatch);
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
                        for (int j = 0; j < sampleIds.length; j++)
                            inputMasks[j][i] = Nd4j.ones(inputShape[0], 1);
                    } else {
                        for (int j = 0; j < sampleIds.length; j++)
                            inputMasks[j][i] = Nd4j.ones(inputShape.clone());
                    }
                }
            }
        }
        if (hasLabelMask) {
            for (int i = 0; i < labelMasks.length; i++) {
                if (labelMasks[i] == null) {
                    int[] labelShape = labels[i].shape();
                    if (labelShape.length == 3) {
                      //  labelMasks[i] = Nd4j.ones(labelShape[0], labelShape[2], 1);
                        //throw new RuntimeException("3D labels should have masks");
                    } else if (labelShape.length == 2 || labelShape.length == 1) {
                        labelMasks[i] = Nd4j.ones(labelShape[0], 1);
                    } else {
                        labelMasks[i] = Nd4j.ones(labelShape.clone());
                    }
                }
            }
        }
        List<MultiDataSet> resultList = new ArrayList<>();
        for (int i = 0; i < sampleIds.length; i++) {
            MultiDataSet result = new org.nd4j.linalg.dataset.MultiDataSet(inputs[i], labels,
                    hasFeatureMask ? inputMasks[i] : null,
                    hasLabelMask ? labelMasks : null);
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
    public void remove() {
        throw new UnsupportedOperationException("Remove is not supported by this iterator.");
    }

    @Override
    public Iterator<List<MultiDataSet>> iterator() {
        reset();
        return this;
    }
}