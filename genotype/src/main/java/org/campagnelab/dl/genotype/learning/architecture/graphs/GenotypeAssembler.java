package org.campagnelab.dl.genotype.learning.architecture.graphs;

import org.campagnelab.dl.framework.domains.DomainDescriptor;
import org.deeplearning4j.nn.conf.ComputationGraphConfiguration;
import org.deeplearning4j.nn.conf.LearningRatePolicy;
import org.deeplearning4j.nn.conf.graph.rnn.DuplicateToTimeSeriesVertex;
import org.deeplearning4j.nn.conf.layers.DenseLayer;
import org.deeplearning4j.nn.conf.layers.GravesLSTM;
import org.deeplearning4j.nn.conf.layers.OutputLayer;
import org.deeplearning4j.nn.conf.layers.RnnOutputLayer;
import org.deeplearning4j.nn.weights.WeightInit;

/**
 * This class provides methods to assemble the isVariant and metaData layers used across genotype assemblers.
 * Created by fac2003 on 1/2/17.
 */
public abstract class GenotypeAssembler {
    protected boolean hasIsVariant;

    protected void appendIsVariantLayer(DomainDescriptor domainDescriptor,
                                        LearningRatePolicy learningRatePolicy,
                                        ComputationGraphConfiguration.GraphBuilder build,
                                        int numIn, WeightInit WEIGHT_INIT,
                                        String lastDenseLayerName) {
        if (hasIsVariant) {
            build.addLayer("isVariant", new OutputLayer.Builder(
                    domainDescriptor.getOutputLoss("isVariant"))
                    .weightInit(WEIGHT_INIT)
                    .activation("softmax").weightInit(WEIGHT_INIT).learningRateDecayPolicy(learningRatePolicy)
                    .nIn(numIn)
                    .nOut(
                            domainDescriptor.getNumOutputs("isVariant")[0]
                    ).build(), lastDenseLayerName);
        }
    }

    protected void appendMetaDataLayer(DomainDescriptor domainDescriptor,
                                       LearningRatePolicy learningRatePolicy,
                                       ComputationGraphConfiguration.GraphBuilder build,
                                       int numIn, WeightInit WEIGHT_INIT,
                                       String lastDenseLayerName) {
        build.addLayer("metaData", new OutputLayer.Builder(
                domainDescriptor.getOutputLoss("metaData"))
                .weightInit(WEIGHT_INIT)
                .activation("softmax").weightInit(WEIGHT_INIT).learningRateDecayPolicy(learningRatePolicy)
                .nIn(numIn)
                .nOut(
                        domainDescriptor.getNumOutputs("metaData")[0]
                ).build(), lastDenseLayerName);
    }

    protected void appendTrueGenotypeLayers(ComputationGraphConfiguration.GraphBuilder build, boolean addTrueGenotypeLabels,
                                            String lastDenseLayerName,
                                            DomainDescriptor domainDescriptor, WeightInit WEIGHT_INIT,
                                            LearningRatePolicy learningRatePolicy,
                                            int numLSTMLayers, int numIn, int numLSTMTrueGenotypeHiddenNodes,
                                            int numLSTMDecoderInputs, int numReductionLayers, float reductionRate) {
        if (addTrueGenotypeLabels) {
            String reductionLayerName = "no layer";
            int numReductionLayerInputNodes = numIn;
            int numReductionLayerOutputNodes = numIn;
            for (int i = 0; i < numReductionLayers; i++) {
                reductionLayerName = "denseReduction_" + i;
                String previousLayerName = i == 0 ? lastDenseLayerName : "denseReduction_" + (i - 1);
                numReductionLayerOutputNodes = (int) (numReductionLayerInputNodes * reductionRate);
                build.addLayer(reductionLayerName, new DenseLayer.Builder()
                        .weightInit(WEIGHT_INIT)
                        .activation("relu")
                        .learningRateDecayPolicy(learningRatePolicy)
                        .nIn(numReductionLayerInputNodes)
                        .nOut(numReductionLayerOutputNodes)
                        .build(), previousLayerName);
                numReductionLayerInputNodes = numReductionLayerOutputNodes;
            }
            build.addVertex("feedForwardLstmDuplicate", new DuplicateToTimeSeriesVertex("trueGenotypeInput"), reductionLayerName);
            String lstmLayerName = "no layer";
            for (int i = 0; i < numLSTMLayers; i++) {
                lstmLayerName = "lstmTrueGenotype_" + i;
                int numLSTMInputNodes = i == 0 ? (numReductionLayerOutputNodes + numLSTMDecoderInputs) : numLSTMTrueGenotypeHiddenNodes;
                String[] layerInputs = i == 0
                        ? new String[]{"trueGenotypeInput", "feedForwardLstmDuplicate"}
                        : new String[]{"lstmTrueGenotype_" + (i - 1)};
                build.addLayer(lstmLayerName, new GravesLSTM.Builder()
                        .activation("softsign")
                        .nIn(numLSTMInputNodes)
                        .nOut(numLSTMTrueGenotypeHiddenNodes)
                        .learningRateDecayPolicy(learningRatePolicy)
                        .weightInit(WEIGHT_INIT)
                        .build(), layerInputs);
            }
            build.addLayer("trueGenotype", new RnnOutputLayer.Builder(domainDescriptor.getOutputLoss("trueGenotype"))
                    .weightInit(WEIGHT_INIT)
                    .activation("softmax")
                    .learningRateDecayPolicy(learningRatePolicy)
                    .nIn(numLSTMTrueGenotypeHiddenNodes)
                    .nOut(domainDescriptor.getNumOutputs("trueGenotype")[0]).build(), lstmLayerName);
        }
    }
}
