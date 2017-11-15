package org.campagnelab.dl.somatic.learning.iterators;

import org.nd4j.linalg.dataset.api.iterator.DataSetIterator;

/**
 * A dataset iterator with a basename.
 * Created by fac2003 on 11/2/16.
 */
public interface NamedDataSetIterator extends DataSetIterator {
    String getBasename();
}
