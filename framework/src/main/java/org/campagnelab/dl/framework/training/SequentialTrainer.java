package org.campagnelab.dl.framework.training;

import it.unimi.dsi.fastutil.ints.IntArrayList;
import it.unimi.dsi.logging.ProgressLogger;
import org.campagnelab.dl.framework.tools.TrainModel;
import org.deeplearning4j.nn.graph.ComputationGraph;
import org.nd4j.linalg.dataset.api.MultiDataSet;
import org.nd4j.linalg.dataset.api.iterator.MultiDataSetIterator;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

/**
 * Traditional sequential training. Works with a single GPU device.
 * Created by fac2003 on 12/1/16.
 */
public class SequentialTrainer implements Trainer {
    private boolean logSpeed;
    static private Logger LOG = LoggerFactory.getLogger(SequentialTrainer.class);

    @Override
    public int train(ComputationGraph computationGraph, MultiDataSetIterator iterator, ProgressLogger progressLogger) {
        int numExamplesUsed = 0;
        int numNanFoundConsecutively = 0;
        while (iterator.hasNext()) {

            MultiDataSet ds = iterator.next();
            // fit the computationGraph:
            computationGraph.fit(ds);

            double score = computationGraph.score();
            if (score != score) {
                // NaN
                numNanFoundConsecutively++;
            } else {
                numNanFoundConsecutively = 0;
            }

            final int numExamples = ds.getFeatures(0).size(0);
            numExamplesUsed += numExamples;

            if (logSpeed) {
                progressLogger.update();
            }
            if (numNanFoundConsecutively > 100) {
                LOG.error("Nan score encountered too many consecutive times");
                return numExamples;
            }
        }
        return numExamplesUsed;
    }

    @Override
    public void setLogSpeed(boolean logSpeed) {
        this.logSpeed = logSpeed;
    }
}
