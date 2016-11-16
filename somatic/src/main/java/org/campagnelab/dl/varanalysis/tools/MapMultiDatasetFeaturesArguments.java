package org.campagnelab.dl.varanalysis.tools;

import org.campagnelab.dl.varanalysis.learning.iterators.MultiDataSetIteratorAdapter;

/**
 * Arguments for MapFeatures.
 * Created by fac2003 on 11/2/16.
 */
public class MapMultiDatasetFeaturesArguments<RecordType> extends MapFeaturesArguments {

    /**
     * Iterators can also be provided programmatically. In this case, they are used rather than creating new ones with
     * the training filenames.
     */
    public MultiDataSetIteratorAdapter adapter =null;


}
