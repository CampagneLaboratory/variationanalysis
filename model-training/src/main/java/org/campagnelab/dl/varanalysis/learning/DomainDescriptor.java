package org.campagnelab.dl.varanalysis.learning;

import org.campagnelab.dl.model.utils.mappers.FeatureMapper;
import org.campagnelab.dl.model.utils.mappers.LabelMapper;
import org.campagnelab.dl.varanalysis.learning.architecture.ComputationalGraphAssembler;
import org.campagnelab.dl.varanalysis.learning.iterators.MultiDataSetRecordIterator;
import org.glassfish.jersey.internal.util.Producer;
import org.nd4j.linalg.lossfunctions.LossFunctions;

import java.util.function.Function;

public abstract class DomainDescriptor<RecordType> {
    public abstract FeatureMapper getFeatureMapper(String inputName);

    public abstract LabelMapper getLabelMapper(String outputName);

    public abstract Function<String[], ? extends MultiDataSetRecordIterator<RecordType>> getIteratorFunction();

    public abstract ComputationalGraphAssembler getComputationalGraph();

    public abstract int[] getNumInputs(String inputName);

    public abstract int[] getNumOutputs(String outputName);

    public abstract int getNumHiddenNodes(String componentName);

    public abstract LossFunctions.LossFunction getOutputLoss(String outputName);

    public int[] getInputShape(int size, String inputName) {
        return getShape(size, () -> getNumInputs(inputName));
    }
    public int[] getLabelShape(int size, String outputName) {
        return getShape(size, () -> getNumOutputs(outputName));
    }

    public int[] getShape(int size, Producer<int[]> p) {
        int[] numInputs =p.call();
        assert numInputs.length <= 2;
        switch (numInputs.length) {
            case 1:
                return new int[]{size, numInputs[0]};
            case 2:
                return new int[]{size, numInputs[0], numInputs[1]};
            default:
                throw new UnsupportedOperationException();
        }
    }
}