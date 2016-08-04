package org.campagnelab.dl.model.utils.mappers;

import org.apache.uima.cas.Feature;
import org.campagnelab.dl.varanalysis.protobuf.BaseInformationRecords;
import org.nd4j.linalg.api.ndarray.INDArray;

/**
 * Created by rct2002 on 7/4/16.
 */
public interface FeatureNameMapper extends FeatureMapper {

    /**
     * Produce the name of a given feature.
     *
     * @param featureIndex The index of the feature to get the name of
     * @return The name of the feature.
     */
    String getFeatureName(int featureIndex);
}
