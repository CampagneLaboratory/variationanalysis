package org.campagnelab.dl.framework.domains.prediction;

import org.campagnelab.dl.framework.domains.DomainDescriptor;
import org.campagnelab.dl.framework.tools.Predict;
import org.campagnelab.dl.framework.tools.VectorReader;
import org.nd4j.linalg.api.ndarray.INDArray;
import org.nd4j.linalg.dataset.api.MultiDataSet;

import java.io.IOException;
import java.util.List;

public class PredictWithVecFile<RecordType> extends PredictWith<RecordType> {

    private final String vecPath;
    private final VectorReader vectorReader;
    private final int numOutputs;

    public PredictWithVecFile(DomainDescriptor domainDescriptor, String vecPath) throws IOException {
        super(domainDescriptor);
        this.vecPath = vecPath;
        String[] outputNames = new String[]{"softmaxGenotype"};
        this.numOutputs = outputNames.length;
        this.vectorReader = new VectorReader(vecPath, 0, outputNames);
    }


    public INDArray[] getModelOutputs(Predict<RecordType> predict, int numOutputs, MultiDataSet dataSet, int batchSize, List<RecordType> records) {
        //assemble outputs from model predictions, and from outputs derived from records (e.g., meta-data):
        INDArray[] outputs = new INDArray[numOutputs];
        INDArray[] modelOutputs = this.vectorReader.getNextBatch(batchSize).getVectors();
        for (int outputIndex = 0; outputIndex < numOutputs; outputIndex++) {
            if (outputIndex < modelOutputs.length) {
                outputs[outputIndex] = modelOutputs[outputIndex];
            } else {
                outputs[outputIndex] = predict.getModelOutput(outputIndex, records);
            }
        }
        return outputs;
    }


}
