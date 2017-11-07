package org.campagnelab.dl.somatic.util;

import it.unimi.dsi.fastutil.longs.Long2BooleanAVLTreeMap;

/**
 * A data structure to remember which genomic sites have already been visited.
 */

public class GenomicSitesVisited {
    Long2BooleanAVLTreeMap visitedSites = new Long2BooleanAVLTreeMap();

    public void visit(int referenceIndex, int position) {
        long value = (((long)(referenceIndex+1))| (((long)(position+1)) << 32));
        visitedSites.put(value, true);

    }

    public boolean wasVisited(int referenceIndex, int position) {
        long value = ((long)(referenceIndex+1)) | ((long)(position+1) << 32);
        return visitedSites.containsKey(value);
    }
}
