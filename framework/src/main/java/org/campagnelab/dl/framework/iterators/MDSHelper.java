package org.campagnelab.dl.framework.iterators;

import org.nd4j.linalg.api.ndarray.INDArray;
import org.nd4j.linalg.dataset.api.MultiDataSet;
import org.nd4j.linalg.factory.Nd4j;

/**
 * Helper for multidatasets. Provides a method to attach an MDS to the computational backend (GPU or CPU).
 * Created by fac2003 on 11/4/17.
 */
public class MDSHelper {

    public static void attach(MultiDataSet ds) {
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
                if (labelsMaskArrays[e] != null) {
                    labelsMaskArrays[e] = labelsMaskArrays[e].detach();
                }
            }
        }
    }
}
