package org.campagnelab.dl.framework.training;

import it.unimi.dsi.logging.ProgressLogger;
import org.deeplearning4j.nn.graph.ComputationGraph;
import org.nd4j.linalg.dataset.api.MultiDataSet;
import org.nd4j.linalg.dataset.api.iterator.MultiDataSetIterator;

/**
 * Traditional sequential training. Works with a single GPU device.
 * Created by fac2003 on 12/1/16.
 */
public class SequentialTrainer implements Trainer {
    @Override
    public int train(ComputationGraph computationGraph, MultiDataSetIterator iterator, ProgressLogger progressLogger) {
        int numExamplesUsed=0;
        while (iterator.hasNext()) {

            MultiDataSet ds = iterator.next();
            // fit the computationGraph:
            computationGraph.fit(ds);
            final int numExamples = ds.getFeatures(0).size(0);
            numExamplesUsed += numExamples;
            progressLogger.lightUpdate();
        }
        return numExamplesUsed;
    }
}
