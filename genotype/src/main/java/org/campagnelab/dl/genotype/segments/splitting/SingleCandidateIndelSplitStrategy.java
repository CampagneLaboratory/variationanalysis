package org.campagnelab.dl.genotype.segments.splitting;

import it.unimi.dsi.fastutil.ints.IntArrayList;
import it.unimi.dsi.fastutil.objects.ObjectArrayList;
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
    private final int windowSize;

    private final int candidateIndelThreshold;

    private int previousCandidateIndel;
    private final BasePositionList basePositionList;
    private boolean verbose;

    public SingleCandidateIndelSplitStrategy(int windowSize, int candidateIndelThreshold, boolean verbose) {
        this.windowSize = windowSize;
        this.candidateIndelThreshold = candidateIndelThreshold;
        this.verbose = verbose;
        basePositionList = new BasePositionList(windowSize);
    }

    /**
     * Applies this strategy to the given argument.
     *
     * @param segment the function argument
     * @return the segment resulting from the splitting.
     */
    @Override
    public List<SingleCandidateIndelSegment> apply(final Segment segment) {
        ObjectArrayList<SingleCandidateIndelSegment> singleCandidateIndelSegments = new ObjectArrayList<>();
        for (BaseInformationRecords.BaseInformation record : segment.getAllRecords()) {
            this.addToBeforeList(record); //keep the before list active
            //complete the sub-segments that are still open
            if (SegmentUtil.hasCandidateIndel(record, this.candidateIndelThreshold) ||
                    SegmentUtil.hasTrueIndel(record)) {
                if ((record.getPosition() - previousCandidateIndel != 1) || (previousCandidateIndel == 0)) { //if they are not consecutive, we open a new subsegment
                    SingleCandidateIndelSegment newSingleCandidateIndelSegment = this.createSubSegment(segment, record);
                    //close the previous subsegments
                    singleCandidateIndelSegments.forEach(SingleCandidateIndelSegment::forceClose);
                    singleCandidateIndelSegments.add(newSingleCandidateIndelSegment);
                } else {
                    //newIndel the previous subsegment end
                    SingleCandidateIndelSegment previous = singleCandidateIndelSegments.get(singleCandidateIndelSegments.size()-1);
                    previous.newIndel(record.getPosition());
                    previous.add(record);
                }
                previousCandidateIndel = record.getPosition();
                
            } else {
                for (SingleCandidateIndelSegment singleCandidateIndelSegment : singleCandidateIndelSegments) {
                    if (singleCandidateIndelSegment.isOpen())
                        singleCandidateIndelSegment.add(record);
                }
            }

        }

        if (verbose && singleCandidateIndelSegments.size() == 0) {
            System.out.println(String.format("No candidate indel found in this segment %d-%d", segment.getFirstPosition(), segment.getLastPosition()));
        }
        //Subsegment as segments
        if (verbose) {
            for (SingleCandidateIndelSegment singleCandidateIndelSegment : singleCandidateIndelSegments) {
                System.out.println(String.format("New subsegment around candidate indel at %d (%d-%d)", singleCandidateIndelSegment.getIndelPosition(),
                        singleCandidateIndelSegment.getFirstPosition(), singleCandidateIndelSegment.getLastPosition()));
            }
            System.out.println(String.format("Subsegments found %d", singleCandidateIndelSegments.size()));
        }

        return singleCandidateIndelSegments;
    }

    private SingleCandidateIndelSegment createSubSegment(Segment parent, BaseInformationRecords.BaseInformation record) {
        return new SingleCandidateIndelSegment(this.basePositionList, parent, record, this.windowSize);
    }

    private void addToBeforeList(BaseInformationRecords.BaseInformation record) {
        this.basePositionList.add(record);
    }

    /**
     * A list that maintains a sized distance between the first and last elements.
     */
    class BasePositionList extends IntArrayList {
        private final long windowSize;
        private int firstElementPosition = 0;
        private int lastElementPosition = 0;

        public BasePositionList(int windowSize) {
            super(windowSize / 2);
            this.windowSize = windowSize;
        }

        public void add(BaseInformationRecords.BaseInformation value) {
            if (value.getPosition() > lastElementPosition || this.size() == 1)
                lastElementPosition = value.getPosition();
            if (value.getPosition() < firstElementPosition || this.size() == 1)
                firstElementPosition = value.getPosition();
            super.add(value.getPosition());
            this.removeEntries();
        }

        /**
         * Removes the entries outside the current window.
         */
        private void removeEntries() {
            for (int i : this) {
                if (lastElementPosition - i > windowSize) {
                    this.rem(i);
                }
            }
        }

        /**
         * {@inheritDoc}
         */
        @Override
        public void clear() {
            super.clear();
            this.firstElementPosition = this.lastElementPosition = 0;
        }
    }
}
