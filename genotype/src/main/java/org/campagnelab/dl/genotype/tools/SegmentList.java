package org.campagnelab.dl.genotype.tools;

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


    protected SegmentList(BaseInformationRecords.BaseInformation from,
                          SequenceSegmentInformationWriter writer, Function<Segment, Segment> function) {
        this.newSegment(from);
        this.function = function;
        this.writer = writer;
    }


    public int getCurrentLocation() {
        return this.currentSegment.getLastPosition();
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

    }

    public String getCurrentLastReferenceId() {
        return this.currentSegment.getLastReferenceId();
    }

    /**
     * Holds the current open segment before it is stored in the list.
     */
    public class Segment {
        private int startPosition = 0;
        private int endPosition = 0;
        List<BaseInformationRecords.BaseInformation> records = new ArrayList<>();

        Segment(BaseInformationRecords.BaseInformation first) {
            //System.out.println("Open a new segment at ref " + first.getReferenceId() + " position " + Integer.toString(first.getPosition()));
            this.records.add(first);
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
        }


        /**
         * Adds a record to the current segment
         *
         * @param record
         */
        protected void add(BaseInformationRecords.BaseInformation record) {
            record.getSamplesList().forEach(sampleInfo -> {
                        /*SegmentInformationRecords.Sample.Builder sampleBuilder = SegmentInformationRecords.Sample.newBuilder();
                        SegmentInformationRecords.Base.Builder baseBuilder = SegmentInformationRecords.Base.newBuilder();
                        //TODO: set real values here
                        baseBuilder.addFeatures(1f);
                        baseBuilder.addLabels(2f);
                        baseBuilder.addTrueLabel("foo");
                        sampleBuilder.addBase(baseBuilder);
                        builder.addSample(sampleBuilder);*/
                    }
            );
            /*SegmentInformationRecords.ReferencePosition.Builder refBuilder = SegmentInformationRecords.ReferencePosition.newBuilder();
            refBuilder.setLocation(record.getPosition());
            refBuilder.setReferenceIndex(record.getReferenceIndex());
            refBuilder.setReferenceId(record.getReferenceId()); */
            this.records.add(record);
        }

        @Override
        public String toString() {
            return "Segment{" +
                    "startPosition=" + this.getFirstPosition() +
                    ", endPosition=" + this.getLastPosition() +
                    ", length=" + actualLength() +
                    '}';
        }

        private int actualLength() {
            if (records.size() == 0)
                return 0;
            else
                return (this.getLastPosition() - this.getFirstPosition()) + 1;
        }

        public int getLastPosition() {
            return records.get(records.size() - 1).getPosition();
        }

        public int getFirstPosition() {
            return records.get(0).getPosition();
        }

        public String getLastReferenceId() {
            return records.get(records.size() - 1).getReferenceId();
        }
    }
}
