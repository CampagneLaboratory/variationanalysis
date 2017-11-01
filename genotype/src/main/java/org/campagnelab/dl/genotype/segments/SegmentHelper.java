package org.campagnelab.dl.genotype.segments;

import org.campagnelab.dl.genotype.segments.splitting.SplitStrategy;
import org.campagnelab.dl.varanalysis.protobuf.BaseInformationRecords;
import org.campagnelab.dl.varanalysis.protobuf.SegmentInformationRecords;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.io.IOException;
import java.nio.charset.Charset;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.*;
import java.util.function.Consumer;
import java.util.function.Function;

import it.unimi.dsi.util.XorShift1024StarRandom;

/**
 * Class to help convert sbi records into segments.
 *
 * @author manuele
 */
public class SegmentHelper {

    private final Function<Segment, Segment> function;
    private final Function<BaseInformationRecords.BaseInformation, SegmentInformationRecords.Base.Builder> fillInFeatures;
    private final SplitStrategy splitStrategy;
    private final Consumer<SegmentInformationRecords.SegmentInformation> segmentConsumer;
    private Segment currentSegment;
    private static Statistics statistics = new Statistics();
    static private Logger LOG = LoggerFactory.getLogger(SegmentHelper.class);
    private boolean collectStatistics;
    private XorShift1024StarRandom random = new XorShift1024StarRandom(28392839);
    private double samplingRate = 0.01;
    int segmentsWithCandidateIndel = 0;
    int segmentsWithTrueIndel = 0;

    /**
     * Creates a new list and a first segment starting from the given record.
     *
     * @param function        the function applied to process the records when the segment is completed.
     * @param fillInFeatures  The function used to fill in features and labels for a post-processed SSI segment.
     * @param segmentConsumer the function called when the segment is completed.
     * @param splitStrategy
     */
    public SegmentHelper(
            Function<Segment, Segment> function,
            Function<BaseInformationRecords.BaseInformation, SegmentInformationRecords.Base.Builder> fillInFeatures,
            Consumer<SegmentInformationRecords.SegmentInformation> segmentConsumer,
            SplitStrategy splitStrategy, boolean collectStatistics) {

        this.function = function;

        this.fillInFeatures = fillInFeatures;
        this.splitStrategy = splitStrategy;
        this.collectStatistics = collectStatistics;
        this.segmentConsumer = segmentConsumer;
    }

    public void setSamplingRate(double samplingRate) {
        this.samplingRate = samplingRate;
    }

    /**
     * Opens a new segment from the record.
     *
     * @param from
     */
    public void newSegment(BaseInformationRecords.BaseInformation from) {
        if (currentSegment != null) {
            this.closeSegment();
            if (collectStatistics) {
                synchronized (statistics) {
                    if (from.getPosition() - currentSegment.getLastPosition() < statistics.minDistance
                            || statistics.minDistance == 0)
                        statistics.minDistance = from.getPosition() - currentSegment.getLastPosition();
                    statistics.addSegment(
                            this.currentSegment.getFirstPosition(), this.currentSegment.getFirstReferenceId(),
                            this.currentSegment.getLastReferenceId(), this.currentSegment.getLastPosition());
                }
            }
        }
        currentSegment = new Segment(fillInFeatures, from);

    }

    /**
     * Closes the current segment.
     */
    private void closeSegment() {
        List<Segment> subSegments = this.splitStrategy.apply(this.currentSegment);

        for (Segment subSegment : subSegments) {
            //sample segments to remove some candidates
            boolean segmentHasCandidateIndel = false;
            boolean segmentHasTrueIndel = false;
            for (BaseInformationRecords.BaseInformation record : subSegment.getAllRecords()) {
                segmentHasCandidateIndel |= SegmentUtil.hasCandidateIndel(record, 0);
                segmentHasTrueIndel |= SegmentUtil.hasTrueIndel(record);

            }

            if (segmentHasCandidateIndel && !segmentHasTrueIndel) {
                // only include some sample of candidates:
                if (random.nextDouble() > samplingRate) {
                    continue;
                }
            }
            segmentsWithCandidateIndel+=segmentHasCandidateIndel?1:0;
            segmentsWithTrueIndel+=segmentHasTrueIndel?1:0;

            //System.out.println(String.format("Processing sub-segment from %d to %d",segment.getFirstPosition(), segment.getLastPosition()));
            try {
                Objects.requireNonNull(this.function);
                Segment processed = this.function.apply(subSegment);
                processed.construct(segmentConsumer);


            } catch (NullPointerException npe) {
                LOG.error("Failed to process segments: ", npe);
                npe.printStackTrace();
                //TODO: throw back this NPE when testing is complete.
                System.out.println(currentSegment);

            } finally {
                this.updateStats();
            }
        }
    }

    public void add(BaseInformationRecords.BaseInformation record) {
        if (currentSegment == null) {
            newSegment(record);
        } else {
            currentSegment.add(record);
        }
    }


    /**
     * Gets the current open segment.
     */
    protected Segment getCurrentSegment() {
        return this.currentSegment;
    }

    /**
     * Close the list.
     */
    public void close() {
        this.closeSegment();
    }

    public void printStats() {
        System.out.printf("Segments with candidate indels: %d%n" +
                "Segments with true indels: %d%n",segmentsWithCandidateIndel,segmentsWithTrueIndel);
        if (!collectStatistics) {
            return;
        }
        System.out.println(statistics);
        Path file = Paths.get(System.currentTimeMillis() + "-segments.tsv");
        List<String> lines = new ArrayList<>();
        List<Long> sortedKeys = new ArrayList<>(statistics.ranges.keySet());
        Collections.sort(sortedKeys);
        int i = 0;
        for (Long start : sortedKeys)
            lines.add(++i + ": " + statistics.ranges.get(start));
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
        if (currentSegment == null) {
            return -1;
        }
        return this.currentSegment.getLastPosition();
    }

    private void updateStats() {
        // TODO we should not synchronize on a global stats. It is quite possible to collect stats in each
        // thread independently and aggregate at the end.
        // However, this is not critical at the moment, so we disable stats collection unless really required.
        if (collectStatistics) {

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

        protected void addSegment(long startIndex, String firstReferenceId, String lastReferenceId, long endIndex) {
            synchronized (ranges) {
                this.ranges.put(startIndex, String.format("%s: %d\t %s:%d", firstReferenceId, startIndex, lastReferenceId, endIndex));
            }
        }
    }


}
