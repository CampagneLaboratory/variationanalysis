package org.campagnelab.dl.genotype.segments;

import org.campagnelab.dl.varanalysis.protobuf.BaseInformationRecords;
import org.campagnelab.dl.varanalysis.protobuf.SegmentInformationRecords;
import org.campagnelab.goby.baseinfo.SequenceSegmentInformationWriter;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.util.Objects;
import java.util.function.Function;

/**
 * Class to help convert sbi records into segments.
 *
 * @author manuele
 */
public class SegmentHelper {

    private final Function<Segment, Segment> function;
    private final Function<BaseInformationRecords.BaseInformation, SegmentInformationRecords.Base.Builder> fillInFeatures;
    private final SequenceSegmentInformationWriter writer;
    private Segment currentSegment;
    private Statistics statistics = new Statistics();
    static private Logger LOG = LoggerFactory.getLogger(SegmentHelper.class);

    /**
     * Creates a new list and a first segment starting from the given record.
     *

     * @param writer         Where the segments will be written when completed.
     * @param function       the function applied to process the records when the segment is completed.
     * @param fillInFeatures The function used to fill in features and labels for a post-processed SSI segment.
     */
    public SegmentHelper(SequenceSegmentInformationWriter writer,
                       Function<Segment, Segment> function,
                       Function<BaseInformationRecords.BaseInformation, SegmentInformationRecords.Base.Builder> fillInFeatures) {

        this.function = function;
        this.writer = writer;
        this.fillInFeatures = fillInFeatures;
    }

    /**
     * Opens a new segment from the record.
     *
     * @param from
     */
    public void newSegment(BaseInformationRecords.BaseInformation from) {
        if (currentSegment != null) {
            this.closeSegment();
            if (from.getPosition() - currentSegment.getLastPosition() < statistics.minDistance
                    || statistics.minDistance == 0)
                statistics.minDistance = from.getPosition() - currentSegment.getLastPosition();
        }
        currentSegment = new Segment(fillInFeatures,from);
        this.add(from);

    }

    /**
     * Closes the current segment.
     */
    private void closeSegment() {
        try {
            Objects.requireNonNull(this.function);
            Segment processed = this.function.apply(currentSegment);
            processed.flush(writer);
            //System.out.println(processed);
        } catch (NullPointerException npe) {
            LOG.error("Failed to process segments: ", npe);
            currentSegment.flush(writer);
            System.out.println(currentSegment);
        } finally {
            this.updateStats();
        }
    }

    public void add(BaseInformationRecords.BaseInformation record) {
        currentSegment.add(record);
    }

    /**
     * Close the list and print the statistics.
     */
    public void close() {
        this.closeSegment();
        this.printStats();
    }

    private void printStats() {
        System.out.println(statistics);
    }

    public int getCurrentReferenceIndex() {
        return this.currentSegment.getLastReferenceIndex();
    }


    public int getCurrentLocation() {
        if (currentSegment==null) {
            return -1;
        }
        return this.currentSegment.getLastPosition();
    }

    private void updateStats() {
        int length = currentSegment.actualLength();
        statistics.totalLength += length;
        if (length > statistics.maxLength)
            statistics.maxLength = length;
        if (length < statistics.minLength || statistics.minLength == 0)
            statistics.minLength = length;
        statistics.numOfSegments++;
    }

    public static boolean isValid(Object record) {
        return false;
    }


    /**
     * Statistics on the list
     */
    class Statistics {
        protected int numOfSegments = 0;
        protected int totalLength = 0;
        protected int minLength = 0;
        protected int maxLength = 0;
        protected int minDistance = 0;

        @Override
        public String toString() {
            return "Statistics{" +
                    "averageLength=" + Math.round(totalLength / numOfSegments) +
                    ", minLength=" + minLength +
                    ", maxLength=" + maxLength +
                    ", minDistance=" + minDistance +
                    '}';
        }
    }



}
