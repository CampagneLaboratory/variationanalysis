package org.campagnelab.dl.somatic.tools;

import org.campagnelab.dl.framework.tools.MapMultiDatasetFeaturesArguments;
import org.campagnelab.dl.somatic.learning.iterators.BaseInformationIterator;

import java.util.List;

/**
 * Arguments for MapFeatures.
 * Created by fac2003 on 11/2/16.
 */
public class SomaticMapFeaturesArguments extends MapMultiDatasetFeaturesArguments {

    /**
     * Iterators can also be provided programmatically. In this case, they are used rather than creating new ones with
     * the training filenames.
     */
    public List<BaseInformationIterator> iterators = null;
}
