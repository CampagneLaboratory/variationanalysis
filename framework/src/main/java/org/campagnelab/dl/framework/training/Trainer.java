package org.campagnelab.dl.framework.training;

import it.unimi.dsi.logging.ProgressLogger;
import org.deeplearning4j.nn.graph.ComputationGraph;
import org.nd4j.linalg.dataset.api.iterator.MultiDataSetIterator;

/**
 * Interface to train computation graphs with an iterator.
 * Created by fac2003 on 12/1/16.
 */
public interface Trainer {
    int train(ComputationGraph graph, MultiDataSetIterator iterator, ProgressLogger pg);
}
