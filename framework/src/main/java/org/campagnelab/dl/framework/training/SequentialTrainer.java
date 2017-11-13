package org.campagnelab.dl.framework.training;

import it.unimi.dsi.fastutil.ints.IntArrayList;
import it.unimi.dsi.logging.ProgressLogger;
import org.campagnelab.dl.framework.iterators.AttachMultiDataSetIterator;
import org.campagnelab.dl.framework.iterators.MDSHelper;
import org.campagnelab.dl.framework.tools.TrainModel;
import org.deeplearning4j.datasets.iterator.AsyncMultiDataSetIterator;
import org.deeplearning4j.nn.conf.WorkspaceMode;
import org.deeplearning4j.nn.graph.ComputationGraph;
import org.nd4j.linalg.api.memory.MemoryWorkspace;
import org.nd4j.linalg.api.memory.conf.WorkspaceConfiguration;
import org.nd4j.linalg.api.memory.enums.AllocationPolicy;
import org.nd4j.linalg.api.memory.enums.LearningPolicy;
import org.nd4j.linalg.api.ndarray.INDArray;
import org.nd4j.linalg.dataset.api.MultiDataSet;
import org.nd4j.linalg.dataset.api.iterator.MultiDataSetIterator;
import org.nd4j.linalg.factory.Nd4j;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

/**
 * Traditional sequential training. Works with a single GPU device.
 * Created by fac2003 on 12/1/16.
 */
public class SequentialTrainer implements Trainer {
    private final WorkspaceConfiguration learningConfig;
    private boolean logSpeed;
    static private Logger LOG = LoggerFactory.getLogger(SequentialTrainer.class);
    double score;
    int n;

    public SequentialTrainer() {
        learningConfig = WorkspaceConfiguration.builder()
                .policyAllocation(AllocationPolicy.STRICT) // <-- this option disables overallocation behavior
                .policyLearning(LearningPolicy.FIRST_LOOP) // <-- this option makes workspace learning after first loop
                .build();
    }


    @Override
    public int train(ComputationGraph computationGraph, MultiDataSetIterator iterator, ProgressLogger progressLogger) {
        int numExamplesUsed = 0;
        int numNanFoundConsecutively = 0;
        score = 0;
        n = 0;


        try (MemoryWorkspace ws = Nd4j.getWorkspaceManager().getAndActivateWorkspace(learningConfig, "TRAINING")) {

            while (iterator.hasNext()) {

                ws.notifyScopeEntered();
                MultiDataSet ds = iterator.next();

                // fit the computationGraph:
                computationGraph.fit(ds);

                double score = computationGraph.score();
                if (score != score) {
                    // NaN
                    numNanFoundConsecutively++;
                } else {
                    numNanFoundConsecutively = 0;
                    this.score += score;
                    this.n++;
                }

                final int numExamples = ds.getFeatures(0).size(0);
                numExamplesUsed += numExamples;
                ds.detach();
                if (logSpeed) {
                    progressLogger.update();
                }
                if (numNanFoundConsecutively > 100) {
                    LOG.error("Nan score encountered too many consecutive times");
                    return numExamples;
                }
                ws.notifyScopeLeft();
            }

        }
        return numExamplesUsed;
    }


    @Override
    public void setLogSpeed(boolean logSpeed) {
        this.logSpeed = logSpeed;
    }

    public double getScore() {
        return score / (double) n;
    }
}

