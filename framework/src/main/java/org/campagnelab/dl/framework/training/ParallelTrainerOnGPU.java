package org.campagnelab.dl.framework.training;

import it.unimi.dsi.logging.ProgressLogger;
import org.deeplearning4j.nn.api.Model;
import org.deeplearning4j.nn.graph.ComputationGraph;
import org.deeplearning4j.optimize.api.IterationListener;
import org.deeplearning4j.optimize.listeners.PerformanceListener;
import org.deeplearning4j.parallelism.ParallelWrapper;
import org.nd4j.linalg.dataset.api.iterator.MultiDataSetIterator;

/**
 * Trainer that trains on multiple GPUs in parallel.
 * Created by fac2003 on 12/1/16.
 */
public class ParallelTrainerOnGPU implements Trainer {
    private final IterationListener perListener = new PerformanceListener(1) {
        private int numNanEncounteredConsecutively;

        @Override
        public void iterationDone(Model model, int iteration) {
            double scoreLocal = model.score();
            if (scoreLocal != scoreLocal) {
                numNanEncounteredConsecutively++;
            } else {
                numNanEncounteredConsecutively = 0;
                score += score;
                n++;
            }
            if (numNanEncounteredConsecutively > 100) {
                wrapper.stopFit();
            }
        }
    };
    ParallelWrapper wrapper;
    int numExamplesPerIterator;
    int miniBatchSize;
    private boolean logSpeed;
    private double score;
    private int n;

    public ParallelTrainerOnGPU(ComputationGraph graph, int miniBatchSize, int totalExamplesPerIterator) {

        String numWorkersString = System.getProperty("framework.parallelWrapper.numWorkers");
        int numWorkers = numWorkersString != null ? Integer.parseInt(numWorkersString) : 4;

        String prefetchBufferString = System.getProperty("framework.parallelWrapper.prefetchBuffer");
        int prefetchBuffer = prefetchBufferString != null ? Integer.parseInt(prefetchBufferString) : 12 * numWorkers;

        String averagingFrequencyString = System.getProperty("framework.parallelWrapper.averagingFrequency");
        int averagingFrequency = averagingFrequencyString != null ? Integer.parseInt(averagingFrequencyString) : 3;

        wrapper = new ParallelWrapper.Builder<>(graph)
                .prefetchBuffer(prefetchBuffer)
                .workers(numWorkers)
                .averagingFrequency(averagingFrequency)
                .reportScoreAfterAveraging(false)
                // .useLegacyAveraging(true)
                .build();
        wrapper.setListeners(perListener);

        this.numExamplesPerIterator = totalExamplesPerIterator;
        this.miniBatchSize = miniBatchSize;
    }

    @Override
    public int train(ComputationGraph graph, MultiDataSetIterator iterator, ProgressLogger pg) {
        score = 0;
        n = 0;
        wrapper.fit(iterator);
        if (logSpeed) {
            pg.update(numExamplesPerIterator);
        }
        return numExamplesPerIterator;
    }

    @Override
    public void setLogSpeed(boolean logSpeed) {
        this.logSpeed = logSpeed;
    }

    @Override
    public double getScore() {
        return score / (double) n;
    }
}
