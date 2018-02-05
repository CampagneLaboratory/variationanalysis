package org.campagnelab.dl.framework.tools;

import org.campagnelab.dl.framework.domains.DomainDescriptor;
import org.campagnelab.dl.framework.domains.prediction.PredictWith;
import org.campagnelab.dl.framework.domains.prediction.PredictionInterpreter;
import org.campagnelab.dl.framework.models.ModelOutputHelper;
import org.deeplearning4j.nn.api.Model;
import org.deeplearning4j.nn.graph.ComputationGraph;
import org.nd4j.linalg.api.ndarray.INDArray;
import org.nd4j.linalg.dataset.api.MultiDataSet;

import java.io.IOException;

/**
 * Helper class to predict with a model and obtain interpreted predictions.
 */
public class PredictWithModel<RecordType> extends PredictWith<RecordType> {
    private final PyTorchModelClient pyTorchModelClient;
    protected DomainDescriptor<RecordType> domainDescriptor;
    ModelOutputHelper outputHelper;
    protected PredictionInterpreter[] interpretors;
    protected Model model;

    public PredictWithModel(DomainDescriptor<RecordType> domainDescriptor, Model model) {
        super(domainDescriptor);
        if (model instanceof PyTorchModel) {
            PyTorchModel pyTorchModel = (PyTorchModel) model;
            PyTorchModelClient.ModelType modelType = pyTorchModel.getPath().startsWith("pytorch:genotyping")
                    ? PyTorchModelClient.ModelType.GENOTYPE
                    : PyTorchModelClient.ModelType.SOMATIC;
            this.pyTorchModelClient = new PyTorchModelClient(pyTorchModel.getLabel(), pyTorchModel.getPath(),
                    modelType, domainDescriptor);
        } else {
            this.pyTorchModelClient = null;
        }
        this.model=model;
    }

    public INDArray[] getModelOutputs(MultiDataSet dataSet, int batchSize) {
        INDArray[] outputPredictions;
        if (pyTorchModelClient != null) {
            try {
                outputPredictions = pyTorchModelClient.predict(dataSet, batchSize);
            } catch (IOException e) {
                throw new RuntimeException("Unable to predict using PyTorch model");
            }
        } else {
            ComputationGraph graph = (ComputationGraph) model;
            outputPredictions = graph.output(false, dataSet.getFeatures());
        }
        return outputPredictions;
    }




}
