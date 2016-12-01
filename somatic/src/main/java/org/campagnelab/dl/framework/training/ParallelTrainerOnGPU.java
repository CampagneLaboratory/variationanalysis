package org.campagnelab.dl.framework.training;

import it.unimi.dsi.logging.ProgressLogger;
import org.deeplearning4j.nn.graph.ComputationGraph;
import org.deeplearning4j.parallelism.ParallelWrapper;
import org.nd4j.linalg.dataset.api.iterator.MultiDataSetIterator;

/**
 * Trainer that trains on multiple GPUs in parallel.
 * Created by fac2003 on 12/1/16.
 */
public class ParallelTrainerOnGPU implements Trainer {
    ParallelWrapper wrapper;
    int numExamplesPerIterator;
    int miniBatchSize;

    public ParallelTrainerOnGPU(ComputationGraph graph, int miniBatchSize, int totalExamplesPerIterator) {

        wrapper = new ParallelWrapper.Builder(graph)
                .prefetchBuffer(12)
                .workers(4)
                .averagingFrequency(1)
                .reportScoreAfterAveraging(false)
                .useLegacyAveraging(false)
                .build();
        this.numExamplesPerIterator = totalExamplesPerIterator;
        this.miniBatchSize = miniBatchSize;
    }

    @Override
    public int train(ComputationGraph graph, MultiDataSetIterator iterator, ProgressLogger pg) {
        System.out.println("Fitting one iterator with paralell wrapper");
        wrapper.fit(iterator);
        pg.update(numExamplesPerIterator);
        return numExamplesPerIterator;
    }
}
