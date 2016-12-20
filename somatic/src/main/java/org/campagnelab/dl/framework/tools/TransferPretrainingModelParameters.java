package org.campagnelab.dl.framework.tools;

import org.campagnelab.dl.framework.architecture.graphs.ComputationGraphAssembler;
import org.campagnelab.dl.framework.domains.DomainDescriptor;
import org.campagnelab.dl.framework.models.ComputationGraphSaver;
import org.campagnelab.dl.framework.models.ModelLoader;
import org.campagnelab.dl.framework.tools.arguments.AbstractTool;
import org.deeplearning4j.nn.api.Model;
import org.deeplearning4j.nn.graph.ComputationGraph;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.io.IOException;
import java.io.UncheckedIOException;

/**
 * Created by joshuacohen on 12/19/16.
 */
public abstract class TransferPretrainingModelParameters<RecordType> extends AbstractTool<TransferPretrainingModelParametersArguments> {
    static private Logger LOG = LoggerFactory.getLogger(TransferPretrainingModelParameters.class);
    private TrainModel<RecordType> trainModel;
    private DomainDescriptor<RecordType> domainDescriptor;
    private ComputationGraphAssembler assembler;
    private ComputationGraph computationGraph;

    @Override
    public TransferPretrainingModelParametersArguments createArguments() {
        return new TransferPretrainingModelParametersArguments();
    }

    @Override
    public void execute() {
        trainModel = createTrainModel();
        domainDescriptor = trainModel.domainDescriptor();
        assembler = domainDescriptor.getComputationalGraph();
        getComputationGraph();
        try {
            transferParams();
        } catch (IOException e) {
            throw new UncheckedIOException(e);
        }
    }

    private void getComputationGraph() {
        ComputationGraphAssembler assembler = domainDescriptor.getComputationalGraph();
        assert assembler != null : "Computational Graph assembler must be defined.";
        assembler.setArguments(trainModel.args());
        for (String inputName : assembler.getInputNames()) {
            int[] domainDescriptorInputs = domainDescriptor.getNumInputs(inputName).clone();
            assert domainDescriptorInputs.length == 2 : "Invalid size for domain descriptor feature";
            boolean padEos = (args().eosIndex != null && args().eosIndex == domainDescriptorInputs[0])
                    || args().eosIndex == null;
            if (padEos) domainDescriptorInputs[0]++;
            assembler.setNumInputs(inputName, domainDescriptorInputs);
        }
        for (String outputName : assembler.getOutputNames()) {
            assembler.setNumOutputs(outputName, domainDescriptor.getNumOutputs(outputName));
            assembler.setLossFunction(outputName, domainDescriptor.getOutputLoss(outputName));
        }
        for (String componentName : assembler.getComponentNames()) {
            assembler.setNumHiddenNodes(componentName, domainDescriptor.getNumHiddenNodes(componentName));
        }

        computationGraph = assembler.createComputationalGraph(domainDescriptor);
        computationGraph.init();
    }

    private void transferParams() throws IOException {
        if (args().pretrainingModelPath != null) {
            ModelLoader pretrainingLoader = new ModelLoader(args().pretrainingModelPath);
            Model savedPretrainingNetwork = pretrainingLoader.loadModel(args().pretrainingModelName);
            ComputationGraph savedPretrainingGraph = savedPretrainingNetwork instanceof ComputationGraph ?
                    (ComputationGraph) savedPretrainingNetwork :
                    null;
            if (savedPretrainingNetwork == null || savedPretrainingGraph == null
                    || savedPretrainingGraph.getUpdater() == null || savedPretrainingGraph.getLayers() == null) {
                LOG.warn("Unable to load model for pretraining from {}", args().pretrainingModelPath);
            } else {
                computationGraph.setUpdater(savedPretrainingGraph.getUpdater());
                for (String inputLayer : assembler.getInputNames()) {
                    computationGraph.getLayer(inputLayer).setParams(
                            savedPretrainingGraph.getLayer(inputLayer).params());
                }
                for (String componentLayer : assembler.getComponentNames()) {
                    computationGraph.getLayer(componentLayer).setParams(
                            savedPretrainingGraph.getLayer(componentLayer).params());
                }
            }
            ComputationGraphSaver graphSaver = new ComputationGraphSaver(args().modelPath);
            graphSaver.saveModel(computationGraph, args().modelPrefix);
        }
    }

    public abstract TrainModel<RecordType> createTrainModel();
}
