package org.campagnelab.dl.varanalysis.learning.iterators;

import org.nd4j.linalg.dataset.DataSet;

/**
 * Created by fac2003 on 7/21/16.
 */
public class SamplingDataset extends org.nd4j.linalg.dataset.DataSet {
    private org.nd4j.linalg.dataset.DataSet delegate;
    private float samplingProbabilities[];

    public SamplingDataset(DataSet delegate, float[] samplingProbabilities, int offsetStartOfMinibatch) {
        this.delegate = delegate;
        this.samplingProbabilities = samplingProbabilities;

        // now
    }
}
