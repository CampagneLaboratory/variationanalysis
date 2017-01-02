package org.campagnelab.dl.genotype.learning.architecture.graphs;

import org.campagnelab.dl.framework.domains.DomainDescriptor;
import org.deeplearning4j.nn.conf.ComputationGraphConfiguration;
import org.deeplearning4j.nn.conf.LearningRatePolicy;
import org.deeplearning4j.nn.conf.layers.OutputLayer;
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
}
