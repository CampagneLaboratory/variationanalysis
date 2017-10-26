package org.campagnelab.dl.genotype.segments.splitting;

import org.campagnelab.dl.genotype.segments.Segment;
import org.campagnelab.dl.genotype.segments.SegmentUtil;
import org.campagnelab.dl.varanalysis.protobuf.BaseInformationRecords;

import java.util.*;

/**
 * Create a sub-segment around each candidate indel.
 */
public class SingleCandidateIndelSplitStrategy implements SplitStrategy {

    /**
     * The window size around the indel.
     */
    private final long windowSize;

    private final int candidateIndelThreshold;

    private final SizedBaseMap beforeIndel;

    public SingleCandidateIndelSplitStrategy(long windowSize, int candidateIndelThreshold) {
        this.windowSize = windowSize;
        this.candidateIndelThreshold = candidateIndelThreshold;
        this.beforeIndel = new SizedBaseMap(windowSize);
    }

    /**
     * Applies this strategy to the given argument.
     *
     * @param segment the function argument
     * @return the segment resulting from the splitting.
     */
    @Override
    public List<Segment> apply(final Segment segment) {
        List<SubSegment> subSegments = new ArrayList<>();

        for (BaseInformationRecords.BaseInformation record : segment.getAllRecords()) {
            this.addToBeforeList(record); //keep the before list active
            if (SegmentUtil.hasCandidateIndel(record, this.candidateIndelThreshold)) {
                SubSegment newSubSegment = this.createSubSegment(segment,record);
                subSegments.add(newSubSegment);
            }
            this.addToBeforeList(record); //keep the before list active for the next subsegments
            //complete the sub-segments that are still open
            for (SubSegment subSegment : subSegments) {
                if (subSegment.isOpen())
                    subSegment.add(record);
            }

        }
         //Subsegment as segments
        List<Segment> segments = new ArrayList<>();
        for (SubSegment subSegment : subSegments){
           segments.add(subSegment.asSegment());
            //System.out.println(String.format("New subsegment around candidate indel at %d (%d-%d)", subSegment.getIndelPosition(),
            //        subSegment.getFirstPosition(), subSegment.getLastPosition()));
        }
        return segments;
    }

    private SubSegment createSubSegment(Segment parent,BaseInformationRecords.BaseInformation record) {
        return new SubSegment(this.beforeIndel,parent,record, this.windowSize);
    }

    private void addToBeforeList(BaseInformationRecords.BaseInformation record) {
        this.beforeIndel.put(record.getPosition(), record);
    }

    /**
     * A map of bases that maintains a sized distance between the first and last elements.
     */
    class SizedBaseMap extends LinkedHashMap<Integer, BaseInformationRecords.BaseInformation> {
        private final long windowSize;
        private int firstElementPosition = 0;
        private int lastElementPosition = 0;

        public SizedBaseMap(long windowSize) {
            super();
            this.windowSize = windowSize;
        }

        /**
         * {@inheritDoc}
         */
        @Override
        public void clear() {
            super.clear();
            this.firstElementPosition = this.lastElementPosition = 0;
        }

        /**
         * Returns <tt>true</tt> if this map should remove its eldest entry.
         *
         * @param eldest The least recently inserted entry in the map
         * @return <tt>true</tt> if the eldest entry should be removed
         * from the map; <tt>false</tt> if it should be retained.
         */
        @Override
        protected boolean removeEldestEntry(Map.Entry<Integer, BaseInformationRecords.BaseInformation> eldest) {
            if (lastElementPosition - eldest.getKey() > windowSize) {
                return this.entrySet().remove(eldest);
            } else {
                return false;
            }

        }

        /**
         * Checks if the list has reached the window size.
         *
         * @return
         */
        protected boolean hasReachedWindowSize() {
            return (lastElementPosition - firstElementPosition >= windowSize);
        }

        /**
         * Associates the specified value with the specified key in this map.
         * If the map previously contained a mapping for the key, the old
         * value is replaced.
         *
         * @param key   key with which the specified value is to be associated
         * @param value value to be associated with the specified key
         * @return the previous value associated with <tt>key</tt>, or
         * <tt>null</tt> if there was no mapping for <tt>key</tt>.
         * (A <tt>null</tt> return can also indicate that the map
         * previously associated <tt>null</tt> with <tt>key</tt>.)
         */
        @Override
        public BaseInformationRecords.BaseInformation put(Integer key, BaseInformationRecords.BaseInformation value) {
            BaseInformationRecords.BaseInformation returnedValue = super.put(key, value);
            if (value.getPosition() > lastElementPosition || this.size() == 1)
                lastElementPosition = value.getPosition();
            if (value.getPosition() < firstElementPosition || this.size() == 1)
                firstElementPosition = value.getPosition();
            return returnedValue;
        }
    }
}
