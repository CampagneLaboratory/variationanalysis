package org.campagnelab.dl.framework.training;

import it.unimi.dsi.fastutil.ints.IntArrayList;
import it.unimi.dsi.logging.ProgressLogger;
import org.campagnelab.dl.framework.tools.TrainModel;
import org.deeplearning4j.datasets.iterator.AsyncMultiDataSetIterator;
import org.deeplearning4j.nn.graph.ComputationGraph;
import org.nd4j.linalg.api.ndarray.INDArray;
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
    double score;
    int n;

    @Override
    public int train(ComputationGraph computationGraph, MultiDataSetIterator iterator, ProgressLogger progressLogger) {
        int numExamplesUsed = 0;
        int numNanFoundConsecutively = 0;
        score = 0;
        n=0;
        //wrap in an async iterator to speed up loading of minibatches to keep the GPU utilized:
        String prefetchBufferString = System.getProperty("framework.parallelWrapper.prefetchBuffer");
        int prefetchBuffer = prefetchBufferString != null ? Integer.parseInt(prefetchBufferString) : 12;
        iterator = new AsyncMultiDataSetIterator(iterator, prefetchBuffer);

        while (iterator.hasNext()) {

            MultiDataSet ds = iterator.next();
            attach(ds);
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

        }
        return numExamplesUsed;
    }

    private void attach(MultiDataSet ds) {
        int e;
        INDArray[] features = ds.getFeatures();
        if (features != null) {
            for (e = 0; e < features.length; ++e) {
                features[e] = features[e].detach();
            }
        }

        INDArray[] labels = ds.getLabels();
        if (labels != null) {
            for (e = 0; e < labels.length; ++e) {
                labels[e] = labels[e].detach();
            }
        }

        INDArray[] featuresMaskArrays = ds.getFeaturesMaskArrays();
        if (featuresMaskArrays != null) {
            for (e = 0; e < featuresMaskArrays.length; ++e) {
                featuresMaskArrays[e] = featuresMaskArrays[e].detach();
            }
        }

        INDArray[] labelsMaskArrays = ds.getLabelsMaskArrays();
        if (labelsMaskArrays != null) {
            for (e = 0; e < labelsMaskArrays.length; ++e) {
                labelsMaskArrays[e] = labelsMaskArrays[e].detach();
            }
        }
    }

    @Override
    public void setLogSpeed(boolean logSpeed) {
        this.logSpeed = logSpeed;
    }

    public double getScore() {
        return score / (double) n;
    }
}
