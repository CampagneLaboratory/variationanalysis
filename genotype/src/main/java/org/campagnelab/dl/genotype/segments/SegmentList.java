package org.campagnelab.dl.genotype.segments;

import org.campagnelab.dl.varanalysis.protobuf.BaseInformationRecords;
import org.campagnelab.dl.varanalysis.protobuf.SegmentInformationRecords;
import org.campagnelab.goby.baseinfo.SequenceSegmentInformationWriter;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.io.IOException;
import java.util.*;
import java.util.function.Function;

/**
 * Segment list for SSI.
 *
 * @author manuele
 */
public class SegmentList {

    private final Function<Segment, Segment> function;
    private final Function<BaseInformationRecords.BaseInformation, SegmentInformationRecords.Base.Builder> fillInFeatures;
    private final SequenceSegmentInformationWriter writer;
    private Segment currentSegment;
    private Statistics statistics = new Statistics();
    static private Logger LOG = LoggerFactory.getLogger(SegmentList.class);

    /**
     * Creates a new list and a first segment starting from the given record.
     *
     * @param from           Initial record in the segment.
     * @param writer         Where the segments will be written when completed.
     * @param function       the function applied to process the records when the segment is completed.
     * @param fillInFeatures The function used to fill in features and labels for a post-processed SSI segment.
     */
    public SegmentList(BaseInformationRecords.BaseInformation from,
                       SequenceSegmentInformationWriter writer,
                       Function<Segment, Segment> function,
                       Function<BaseInformationRecords.BaseInformation, SegmentInformationRecords.Base.Builder> fillInFeatures) {
        this.newSegment(from);
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
        currentSegment = new Segment(from);
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

    /**
     * Holds the current open segment before it is stored in the list.
     */
    public class Segment {
        private int firstPosition = 0;
        private int firstReferenceIndex = 0;
        private String firstReferenceId = "";
        private int lastPosition = 0;
        private String lastReferenceId = "";
        private int lastReferenceIndex = 0;
        public RecordList recordList = new RecordList();

        Segment(BaseInformationRecords.BaseInformation first) {
            //System.out.println("Open a new segment at ref " + first.getReferenceId() + " position " + Integer.toString(first.getPosition()));
            this.add(first);
            this.firstPosition = first.getPosition();
            this.firstReferenceIndex = first.getReferenceIndex();
            this.firstReferenceId = first.getReferenceId();
        }

        private void setAsLast(BaseInformationRecords.BaseInformation record) {
            this.lastPosition = record.getPosition();
            this.lastReferenceId = record.getReferenceId();
            this.lastReferenceIndex = record.getReferenceIndex();
        }

        /**
         * Flushes the segment in the SBI writer.
         *
         * @param writer
         */
        protected void flush(SequenceSegmentInformationWriter writer) {

            SegmentInformationRecords.SegmentInformation.Builder builder = SegmentInformationRecords.SegmentInformation.newBuilder();
            SegmentInformationRecords.ReferencePosition.Builder refBuilder = SegmentInformationRecords.ReferencePosition.newBuilder();
            refBuilder.setLocation(this.getFirstPosition());
            refBuilder.setReferenceIndex(this.getFirstReferenceIndex());
            refBuilder.setReferenceId(this.getFirstReferenceId());
            builder.setStartPosition(refBuilder.build());
            refBuilder = SegmentInformationRecords.ReferencePosition.newBuilder();
            refBuilder.setLocation(this.getLastPosition());
            refBuilder.setReferenceIndex(this.getLastReferenceIndex());
            refBuilder.setReferenceId(this.getLastReferenceId());
            builder.setEndPosition(refBuilder.build());
            builder.setLength(actualLength());
            final long[] segmentStats = {0L, 0L, 0L};
            recordList.forEach(record -> {
                record.getSamplesList().forEach(sample -> {
                    SegmentInformationRecords.Sample.Builder sampleBuilder = SegmentInformationRecords.Sample.newBuilder();

                    SegmentInformationRecords.Base.Builder base = fillInFeatures.apply(record);
                    sampleBuilder.addBase(base);
                    builder.addSample(sampleBuilder);
                    segmentStats[0]++;
                    //System.out.println("New base " + segmentStats[0] );
                    segmentStats[1] = base.getFeaturesCount();
                    segmentStats[2] = base.getLabelsCount();
                });
            });
            try {
                writer.appendEntry(builder.build());
                writer.setEntryBases(segmentStats[0]);
                writer.setEntryLabels(segmentStats[1]);
                writer.setEntryFeatures(segmentStats[2]);
            } catch (IOException e) {
                e.printStackTrace();
            }
        }

        /**
         * Adds a record to the current segment
         *
         * @param record
         */
        protected void add(BaseInformationRecords.BaseInformation record) {
            this.recordList.add(record);
            this.setAsLast(record);
        }

        @Override
        public String toString() {
            return String.format("Segment{start=%s:%d end=%s:%d length=%d}%n", getFirstReferenceId(), getFirstPosition(),
                    getLastReferenceId(), getLastPosition(), actualLength());
        }

        private int actualLength() {
            if (recordList.size() == 0)
                return 0;
            else
                return (this.getLastPosition() - this.getFirstPosition()) + 1;
        }

        public String getFirstReferenceId() {
            return this.firstReferenceId;
        }

        public int getFirstPosition() {
            return this.firstPosition;
        }

        public int getFirstReferenceIndex() {
            return this.firstReferenceIndex;
        }


        public String getLastReferenceId() {
            return this.lastReferenceId;
        }

        public int getLastPosition() {
            return this.lastPosition;
        }

        public int getLastReferenceIndex() {
            return this.lastReferenceIndex;
        }

        public boolean hasTrueGenotype(String trueGenotype) {
            for (BaseInformationRecords.BaseInformation record : recordList) {
                if (record.getTrueGenotype().equals(trueGenotype)) {
                    return true;
                }

            }
            return false;
        }


    }

}
