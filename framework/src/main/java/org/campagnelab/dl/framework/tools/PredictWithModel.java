package org.campagnelab.dl.framework.tools;

import org.campagnelab.dl.framework.domains.DomainDescriptor;
import org.campagnelab.dl.framework.domains.prediction.PredictWith;
import org.campagnelab.dl.framework.domains.prediction.PredictionInterpreter;
import org.campagnelab.dl.framework.models.ModelOutputHelper;
import org.deeplearning4j.nn.api.Model;
import org.deeplearning4j.nn.graph.ComputationGraph;
import org.nd4j.linalg.api.ndarray.INDArray;
import org.nd4j.linalg.dataset.api.MultiDataSet;

/**
 * Helper class to predict with a model and obtain interpreted predictions.
 */
public class PredictWithModel<RecordType> extends PredictWith<RecordType> {
    protected DomainDescriptor<RecordType> domainDescriptor;
    ModelOutputHelper outputHelper;
    protected PredictionInterpreter[] interpretors;
    protected Model model;

    public PredictWithModel(DomainDescriptor<RecordType> domainDescriptor, Model model) {
        super(domainDescriptor);
        this.model=model;
    }

    public INDArray[] getModelOutputs(MultiDataSet dataSet, int batchSize) {
        ComputationGraph graph=(ComputationGraph)model;
        INDArray[] outputPredictions = graph.output(false,dataSet.getFeatures());
        return outputPredictions;
    }




}
