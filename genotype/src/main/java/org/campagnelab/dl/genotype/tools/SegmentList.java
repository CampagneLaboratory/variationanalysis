package org.campagnelab.dl.genotype.tools;

import org.campagnelab.dl.varanalysis.protobuf.BaseInformationRecords;

import java.util.Iterator;

/**
 * Segment list for SSI.
 *
 * @author manuele
 */
public class SegmentList implements Iterable<SegmentList.Segment>{

    private Segment currentSegment;

    private int currentLastPosition = 0;
    private int currentLastReferenceIndex = 0;
    private String currentLastReferenceId = "";


    protected SegmentList(BaseInformationRecords.BaseInformation from) {
        this.newSegment(from);
    }


    public int getCurrentLocation() {
        return currentLastPosition;
    }

    /**
     * Opens a new segment from the record.
     * @param from
     */
    public void newSegment(BaseInformationRecords.BaseInformation from) {
        if (currentSegment != null)
            this.closeSegment();
        currentSegment = new Segment(from);
        this.add(from);

    }

    public void closeSegment() {
       currentSegment.close();
       currentSegment.printStats();
    }

    public void add(BaseInformationRecords.BaseInformation record) {
       currentSegment.add(record);
       setAsLast(record);
    }

    /**
     * Sets the record as the last one in the current segment.
     * @param record
     */
    private void setAsLast(BaseInformationRecords.BaseInformation record) {
        currentLastPosition = record.getPosition();
        currentLastReferenceIndex = record.getReferenceIndex();
        currentLastReferenceId = record.getReferenceId();
    }

    /**
     * Returns an iterator over elements of type {@code T}.
     *
     * @return an Iterator.
     */
    @Override
    public Iterator<Segment> iterator() {
        return null;
    }

    public void close() {
        this.closeSegment();
        this.printStats();
    }

    private void printStats() {
        
    }

    /**
     * Holds the current open segment before it is stored in the list.
     */
    class Segment {
        private int startPosition = 0;
        private int endPosition = 0;

        Segment(BaseInformationRecords.BaseInformation first) {
            System.out.println("Open a new segment at position " + Integer.toString(first.getPosition()));
            this.startPosition = first.getPosition();
            this.endPosition = first.getPosition();
        }

        protected void close() {
            //if (builder != null) {
                //close the previous segment
                System.out.println("Close the segment.");
                 /*try {
                   writer.appendEntry(builder.build());
                    //set the current* as end position
                    SegmentInformationRecords.ReferencePosition.Builder refBuilder = SegmentInformationRecords.ReferencePosition.newBuilder();
                    refBuilder.setLocation(currentLastPosition);
                    refBuilder.setReferenceIndex(currentLastReferenceIndex);
                    refBuilder.setReferenceId(currentLastReferenceId);
                    builder.setEndPosition(refBuilder);
                    printSegmentStats();
                } catch (IOException e) {
                    System.err.println("Unable to close the previous segment");
                    e.printStackTrace();
                } finally {
                    builder = null;
                    currentLastPosition = 0;
                    currentLastReferenceId = "";
                    currentLastReferenceIndex = 0;
                }   */
                //create statistics here.
            //}
        }

        protected void printStats() {
            System.out.println("Start position:" + startPosition);
            System.out.println("End position:" + endPosition);
            System.out.println("Length:" + (endPosition - startPosition + 1));
        }

        /**
         * Adds a record to the current segment
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
            this.endPosition = record.getPosition();

        }
    }
}
