package org.campagnelab.dl.framework.domains.prediction;

import org.campagnelab.dl.framework.domains.DomainDescriptor;
import org.campagnelab.dl.framework.tools.VectorReader;
import org.campagnelab.dl.framework.tools.VectorWriter;
import org.nd4j.linalg.api.ndarray.INDArray;
import org.nd4j.linalg.dataset.api.MultiDataSet;
import org.nd4j.linalg.factory.Nd4j;

import java.io.IOException;

public class PredictWithVecFile<RecordType> extends PredictWith<RecordType> {

    private final String vecPath;
    private final VectorReader vectorReader;
    private final int numOutputs;

    public PredictWithVecFile(DomainDescriptor domainDescriptor, String vecPath) throws IOException {
        super(domainDescriptor);
        this.vecPath = vecPath;
        String[] outputNames = {"softmaxGenotype","metadata"};// TODO: obtain from domainDescriptor.getOutputNames();
        this.numOutputs = outputNames.length;
        this.vectorReader = new VectorReader(vecPath, 0, outputNames);
    }


    public INDArray[] getModelOutputs(MultiDataSet dataSet, int batchSize) {

        return this.vectorReader.getNextBatch(batchSize).getVectors();
        /*INDArray[] outputPredictions = new INDArray[numOutputs];
        int[] shape;

        final VectorWriter.VectorProperties.VectorPropertiesVector[] vectors = this.vectorReader.getVectorProperties().getVectors();
        for (int exampleIndex = 0; exampleIndex < batchSize; exampleIndex++) {
            VectorReader.RecordVectors outputs = this.vectorReader.getNextExample();

            for (int outputIndex = 0; outputIndex < numOutputs; outputIndex++) {
                if (outputPredictions[outputIndex]==null) {

                    outputPredictions[outputIndex]= Nd4j.createUninitializedDetached(vectors[outputIndex].getVectorDimension())
                }
                outputPredictions[outputIndex].getRow(exampleIndex).assign(outputs.getVectors()[outputIndex]);
            }
        }
        return outputPredictions;*/
    }


}
