package org.campagnelab.dl.genotype.tools;

import it.unimi.dsi.fastutil.objects.ObjectArrayList;
import org.campagnelab.dl.varanalysis.protobuf.BaseInformationRecords;
import org.campagnelab.dl.varanalysis.protobuf.SegmentInformationRecords;
import org.campagnelab.goby.baseinfo.SequenceSegmentInformationWriter;

import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.function.Function;

/**
 * Segment list for SSI.
 *
 * @author manuele
 */
public class SegmentList {

    private final Function<Segment, Segment> function;
    private final SequenceSegmentInformationWriter writer;
    private Segment currentSegment;
    private Statistics statistics = new Statistics();

    /**
     * Creates a new list and a first segment starting from the given record.
     * @param from
     * @param writer
     * @param function the function applied to the records when the segment is completed.
     */
    protected SegmentList(BaseInformationRecords.BaseInformation from,
                          SequenceSegmentInformationWriter writer, Function<Segment, Segment> function) {
        this.newSegment(from);
        this.function = function;
        this.writer = writer;
    }

    /**
     * Opens a new segment from the record.
     *
     * @param from
     */
    public void newSegment(BaseInformationRecords.BaseInformation from) {
        if (currentSegment != null)
            this.closeSegment();
        currentSegment = new Segment(from);
        this.add(from);

    }

    public void closeSegment() {
        if (this.function != null) {
            Segment processed = this.function.apply(currentSegment);
            processed.flush(writer);
            System.out.println(processed);
        } else {
            currentSegment.flush(writer);
            System.out.println(currentSegment);
        }
    }

    public void add(BaseInformationRecords.BaseInformation record) {
        currentSegment.add(record);
    }

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

    /**
     * Statistics on the list
     */
    class Statistics {
        protected int numOfSegments = 0;
        protected int totalLength = 0;
        protected int minLength = 0;
        protected int maxLength = 0;

        @Override
        public String toString() {
            return "Statistics{" +
                    "averageLength=" + Math.round(totalLength/numOfSegments) +
                    ", minLength=" + minLength +
                    ", maxLength=" + maxLength +
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
        private String lastReferenceId= "";
        private int lastReferenceIndex = 0;
        List<BaseInformationRecords.BaseInformation> records = new ObjectArrayList<>();

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
            refBuilder.setLocation(records.get(0).getPosition());
            refBuilder.setReferenceIndex(records.get(0).getReferenceIndex());
            refBuilder.setReferenceId(records.get(0).getReferenceId());
            builder.setStartPosition(refBuilder.build());
            refBuilder = SegmentInformationRecords.ReferencePosition.newBuilder();
            refBuilder.setLocation(records.get(records.size() - 1).getPosition());
            refBuilder.setReferenceIndex(records.get(records.size() - 1).getReferenceIndex());
            refBuilder.setReferenceId(records.get(records.size() - 1).getReferenceId());
            builder.setEndPosition(refBuilder.build());
            builder.setLength(actualLength());
            records.forEach(record -> {
                    record.getSamplesList().forEach(sample -> {
                        SegmentInformationRecords.Sample.Builder sampleBuilder = SegmentInformationRecords.Sample.newBuilder();
                        SegmentInformationRecords.Base.Builder baseBuilder = SegmentInformationRecords.Base.newBuilder();
                        //TODO: set real values here
                        baseBuilder.addFeatures(1f);
                        baseBuilder.addLabels(2f);
                        baseBuilder.addTrueLabel("foo");
                        sampleBuilder.addBase(baseBuilder);
                        builder.addSample(sampleBuilder);
                    });
            });
            try {
                writer.appendEntry(builder.build());
            } catch (IOException e) {
                e.printStackTrace();
            }
            int length = this.actualLength();
            statistics.totalLength += length;
            if (length > statistics.maxLength)
                statistics.maxLength = length;
            if (length < statistics.minLength || statistics.minLength == 0)
                statistics.minLength = length;
            statistics.numOfSegments++;
        }


        /**
         * Adds a record to the current segment
         *
         * @param record
         */
        protected void add(BaseInformationRecords.BaseInformation record) {
            this.records.add(record);
            this.setAsLast(record);
        }

        @Override
        public String toString() {
            return String.format("Segment{start=%s:%d end=%s:%d length=%d}%n", getFirstReferenceId(), getFirstPosition(),
                    getLastReferenceId(), getLastPosition(),actualLength());
        }

        private int actualLength() {
            if (records.size() == 0)
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

        public String getLastReferenceId() {
            return this.lastReferenceId;
        }
        public int getLastPosition() {
            return this.lastPosition;
        }

        public int getLastReferenceIndex() {
            return this.lastReferenceIndex;
        }
    }
}
