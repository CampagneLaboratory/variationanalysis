package org.campagnelab.dl.somatic.mappers;

import org.campagnelab.dl.framework.mappers.ConcatFeatureMapper;
import org.campagnelab.dl.framework.mappers.FeatureMapper;
import org.campagnelab.dl.framework.mappers.FeatureNameMapper;

import java.util.Arrays;

/**
 * Created by rct66 on 8/4/16.
 *
 * This class is an alternative to ConcatFeatureMapper which enables featurename extraction.
 * This was created so that deprecated mappers would not need to have feature naming implemented.
 */
public class NamingConcatFeatureMapper<RecordType> extends ConcatFeatureMapper<RecordType> implements FeatureNameMapper<RecordType> {


    @SafeVarargs
    public NamingConcatFeatureMapper(FeatureNameMapper<RecordType>... featureMappers) {
        super(featureMappers);
    }

    public String getFeatureName(int i) {
        int indexOfDelegate = Arrays.binarySearch(offsets, i);
        if (indexOfDelegate < 0) {
            indexOfDelegate = -(indexOfDelegate + 1) - 1;
        }
        return ((FeatureNameMapper)this.mappers[indexOfDelegate]).getFeatureName(i - offsets[indexOfDelegate]);
    }


}
