package org.campagnelab.dl.genotype.tools;

import it.unimi.dsi.fastutil.objects.ObjectArrayList;
import org.apache.commons.lang3.tuple.MutablePair;
import org.campagnelab.dl.varanalysis.protobuf.SegmentInformationRecords;

/**
 * List of bases accumulated when analysing a prediction on segments.
 *
 * @author manuele
 */
public class VCFLine extends ObjectArrayList<VCFLine.IndexedBase> {

    private boolean isIndel = false;
    private int lastBaseLocation = 0;
    private int lastGapLocation = 0;
    private boolean noAnchorBeforeGap = false;


    /**
     * Detects according to the next base, if the line needs to be written in the VCF
     *
     * @param nextBase the candidate next base
     * @return
     */
    public boolean needToFlush(IndexedBase nextBase) {
        if (this.isEmpty())
            return false;
        if (!isIndel)
            //check if they are contiguous
            return lastBaseLocation != nextBase.getKey().getLocation();
        else {
            //the next base after a gap is always accepted
            if (lastBaseLocation == lastGapLocation)
                return false;
            else {
                //check if they are contiguous
                return lastBaseLocation != nextBase.getKey().getLocation();
            }
        }

    }

    /**
     * Marks the current content of the queue as an indel
     */
    private void markAsIndel(IndexedBase base) {
        this.isIndel = true;
        this.lastGapLocation = base.getKey().getLocation();
    }

    @Override
    public void add(int index, IndexedBase base) {
        super.add(index, base);
        this.lastBaseLocation = base.getKey().getLocation();
    }

    @Override
    public boolean add(IndexedBase base) {
        this.lastBaseLocation = base.getKey().getLocation();
        return super.add(base);
    }

    @Override
    public void clear() {
        super.clear();
        this.isIndel = false;
        this.lastBaseLocation = 0;
        this.lastGapLocation = 0;
    }

    public boolean isIndel() {
        return isIndel;
    }

    public void addBaseWithPredictedGap(IndexedBase indexedBase) {
        if (this.isEmpty()) {
            noAnchorBeforeGap = true;
        }
        this.add(indexedBase);
        this.markAsIndel(indexedBase);
    }


    public boolean hasNoAnchorBeforeGap() {
        return noAnchorBeforeGap;
    }


    /**
     * Base with its current position in the segment.
     */
    public static class IndexedBase extends MutablePair<SegmentInformationRecords.Base, Integer> {

        public IndexedBase(SegmentInformationRecords.Base base, Integer index) {
            super(base, index);
        }
    }


}
