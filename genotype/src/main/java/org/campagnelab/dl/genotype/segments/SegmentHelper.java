package org.campagnelab.dl.genotype.segments;

import org.campagnelab.dl.varanalysis.protobuf.BaseInformationRecords;
import org.campagnelab.dl.varanalysis.protobuf.SegmentInformationRecords;
import org.campagnelab.goby.baseinfo.SequenceSegmentInformationWriter;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.io.IOException;
import java.nio.charset.Charset;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.*;
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
    private static  Statistics statistics = new Statistics();
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
            synchronized (statistics) {
                if (from.getPosition() - currentSegment.getLastPosition() < statistics.minDistance
                        || statistics.minDistance == 0)
                    statistics.minDistance = from.getPosition() - currentSegment.getLastPosition();
                statistics.addSegment(this.currentSegment.getFirstPosition(), this.currentSegment.getLastPosition());
            }
        }
        currentSegment = new Segment(fillInFeatures,from);

    }

    /**
     * Closes the current segment.
     */
    private void closeSegment() {
        try {
            Objects.requireNonNull(this.function);
            Segment processed = this.function.apply(currentSegment);
            processed.flush(writer);

            if (processed.actualLength()>300) {
                System.out.println(processed);
                System.out.println("STOP");
            }
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
     * Close the list.
     */
    public void close() {
        this.closeSegment();
    }

    public void printStats() {

        System.out.println(statistics);
        Path file = Paths.get(System.currentTimeMillis() + "-segments.tsv");
        List<String> lines = new ArrayList<>();
        List<Long> sortedKeys=new ArrayList<>(statistics.ranges.keySet());
        Collections.sort(sortedKeys);
        int i = 0;
        for (Long start : sortedKeys )
            lines.add(++i + ": " + statistics.ranges.get(start)) ;
        try {
            Files.write(file, lines, Charset.forName("UTF-8"));
        } catch (IOException e) {
            e.printStackTrace();
        }
        //statistics.ranges.clear();
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
        synchronized (statistics) {
            int length = currentSegment.actualLength();
            statistics.totalLength += length;
            if (length > statistics.maxLength)
                statistics.maxLength = length;
            if (length < statistics.minLength || statistics.minLength == 0)
                statistics.minLength = length;
            statistics.numOfSegments++;
        }
    }

    public static boolean isValid(Object record) {
        return false;
    }


    /**
     * Statistics on the list
     */
    static class Statistics {
        protected int numOfSegments = 0;
        protected int totalLength = 0;
        protected int minLength = 0;
        protected int maxLength = 0;
        protected int minDistance = 0;
        protected static Map<Long, String> ranges = new HashMap<>();

        @Override
        public String toString() {
            return "Statistics{" +
                    "averageLength=" + Math.round(totalLength / numOfSegments) +
                    ", minLength=" + minLength +
                    ", maxLength=" + maxLength +
                    ", minDistance=" + minDistance +
                    '}';
        }

        protected void addSegment(long startIndex, long endIndex) {
            synchronized (ranges) {
                this.ranges.put(startIndex, String.format("%d\t%d", startIndex, endIndex));
            }
        }
    }



}
