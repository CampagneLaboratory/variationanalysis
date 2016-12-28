package org.campagnelab.dl.framework.tools;

import org.apache.commons.io.FileUtils;
import org.campagnelab.dl.framework.architecture.graphs.ComputationGraphAssembler;
import org.campagnelab.dl.framework.domains.DomainDescriptor;
import org.campagnelab.dl.framework.models.ComputationGraphSaver;
import org.campagnelab.dl.framework.models.ModelLoader;
import org.campagnelab.dl.framework.tools.arguments.AbstractTool;
import org.deeplearning4j.nn.api.Model;
import org.deeplearning4j.nn.graph.ComputationGraph;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.io.*;
import java.util.Date;
import java.util.Properties;

/**
 * Created by joshuacohen on 12/19/16.
 */
public abstract class TransferPretrainingModelParameters<RecordType> extends AbstractTool<TransferPretrainingModelParametersArguments> {
    static private Logger LOG = LoggerFactory.getLogger(TransferPretrainingModelParameters.class);
    private DomainDescriptor<RecordType> domainDescriptor;
    private ComputationGraphAssembler assembler;
    private ComputationGraph computationGraph;
    private String modelPath;

    @Override
    public void execute() {
        try {
            modelPath = args().modelPath != null ? args().modelPath : "models/" + Long.toString(new Date().getTime());
            FileUtils.forceMkdir(new File(modelPath));
            domainDescriptor = (DomainDescriptor<RecordType>) Class.forName(args().domainDescriptorName())
                    .getConstructor(String.class)
                    .newInstance(args().pretrainingModelPath);
            assembler = domainDescriptor.getComputationalGraph();
            getComputationGraph();
            transferParams();
            transferProperties();
        } catch (Exception e) {
            throw new RuntimeException("Couldn't transfer parameters", e);
        }
    }

    private void getComputationGraph() {
        ComputationGraphAssembler assembler = domainDescriptor.getComputationalGraph();
        assert assembler != null : "Computational Graph assembler must be defined.";
        assembler.setArguments(args());
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

    private void transferProperties() throws IOException {
        Properties configProperties = new Properties();
        Properties domainProperties = new Properties();
        Reader configReader = new FileReader(new File(args().pretrainingModelPath, "config.properties"));
        Reader domainReader = new FileReader(new File(args().pretrainingModelPath, "domain.properties"));
        configProperties.load(configReader);
        domainProperties.load(domainReader);
        Writer configWriter = new FileWriter(new File(modelPath, "config.properties"));
        Writer domainWriter = new FileWriter(new File(modelPath, "domain.properties"));
        configProperties.store(configWriter,
                String.format("Config properties created via pretraining parameter transfer from %s",
                        args().pretrainingModelPath));
        domainProperties.store(domainWriter,
                String.format("Domain properties created via pretraining parameter transfer from %s",
                        args().pretrainingModelPath));
        configWriter.close();
        domainWriter.close();
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
                throw new RuntimeException(String.format("Unable to load model for pretraining from %s",
                        args().pretrainingModelPath));
            } else {
                for (String inputLayer : assembler.getInputNames()) {
                    computationGraph.getLayer(inputLayer).setParams(
                            savedPretrainingGraph.getLayer(inputLayer).params());
                }
                for (String componentLayer : assembler.getComponentNames()) {
                    computationGraph.getLayer(componentLayer).setParams(
                            savedPretrainingGraph.getLayer(componentLayer).params());
                }
            }
            String modelPrefix = args().modelPrefix != null ? args().modelPrefix : "pretraining";
            ComputationGraphSaver graphSaver = new ComputationGraphSaver(modelPath);
            graphSaver.saveModel(computationGraph, modelPrefix);
        }
    }
}
