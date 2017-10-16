package org.campagnelab.dl.genotype.tools;

import it.unimi.dsi.fastutil.floats.FloatArrayList;
import it.unimi.dsi.fastutil.ints.Int2ObjectAVLTreeMap;
import it.unimi.dsi.fastutil.ints.Int2ObjectMap;
import it.unimi.dsi.fastutil.objects.Object2ObjectArrayMap;
import it.unimi.dsi.fastutil.objects.ObjectArrayList;
import org.campagnelab.dl.varanalysis.protobuf.BaseInformationRecords;
import org.campagnelab.dl.varanalysis.protobuf.SegmentInformationRecords;
import org.campagnelab.goby.baseinfo.SequenceSegmentInformationWriter;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import scala.Char;

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
     * @param from
     * @param writer
     * @param function the function applied to the records when the segment is completed.
     */
    protected SegmentList(BaseInformationRecords.BaseInformation from,
                          SequenceSegmentInformationWriter writer, Function<Segment, Segment> function,
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
            final long[] segmentStats = {0L, 0L ,0L};
            recordList.forEach(record -> {
                    record.getSamplesList().forEach(sample -> {
                        SegmentInformationRecords.Sample.Builder sampleBuilder = SegmentInformationRecords.Sample.newBuilder();

                        SegmentInformationRecords.Base.Builder base =fillInFeatures.apply(record);
                        sampleBuilder.addBase(base);builder.addSample(sampleBuilder);
                    segmentStats[0]++;
                        //System.out.println("New base " + segmentStats[0] );
                        segmentStats[1] = base.getFeaturesCount();
                        segmentStats[2] = base.getLabelsCount();});
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


        class GenotypesAtSite {
            //one true genotype per channel (#channels==ploidy)
            char[] trueGenotypes;

        }

        class RecordList implements Iterable<BaseInformationRecords.BaseInformation> {

            ObjectArrayList<BaseInformationRecords.BaseInformation> records = new ObjectArrayList<>();

            /**
             * Returns an iterator over elements of type {@code T}.
             *
             * @return an Iterator.
             */
            @Override
            public Iterator<BaseInformationRecords.BaseInformation> iterator() {
                return records.iterator();
            }


            public int size() {
                return records.size();
            }

            public void add(BaseInformationRecords.BaseInformation record) {
                records.add(record);
            }

            public BaseInformationRecords.BaseInformation insertAfter(BaseInformationRecords.BaseInformation previous,
                                                                      BaseInformationRecords.BaseInformation buildFrom,
                                                                      char insertedDeleted, int offset) {
                BaseInformationRecords.BaseInformation.Builder copy = buildFrom.toBuilder();

                String trueFrom = previous.getTrueFrom();
                if (trueFrom.length() > offset + 1) {
                    copy.setTrueFrom(previous.getTrueFrom().substring(offset, offset + 1));
                } else {
                    copy.setTrueFrom("");
                }
                copy.setTrueGenotype(Character.toString(insertedDeleted));
                copy = adjustCounts(copy, offset);
                final BaseInformationRecords.BaseInformation builtCopy = copy.build();

                records.add(records.indexOf(previous) + 1, builtCopy);
                return previous;

            }

            /**
             * Adjust counts to reduce indels to a single base, using offset and the current indel sequences to determine where
             * to increase the base count.
             *
             * @param copy   record that needs to have counts adjusted.
             * @param offset offset inside the indel sequence, which identifies the base to increment.
             * @return a builder where the adjustment has been made.
             */
            private BaseInformationRecords.BaseInformation.Builder adjustCounts(BaseInformationRecords.BaseInformation.Builder copy, int offset) {
                // we store counts in a map for easy access (map keyed on to sequence of the count):
                Map<String, BaseInformationRecords.CountInfo.Builder> counts = new Object2ObjectArrayMap<>();
                // we know we may need a gap count, so we add one, because none in the sbi:

                int sampleIndex = 0;
                int countIndex = 0;
                boolean needsGap = true;
                for (BaseInformationRecords.SampleInfo.Builder sample : copy.getSamplesBuilderList()) {
                    for (BaseInformationRecords.CountInfo.Builder count : sample.getCountsBuilderList()) {
                        String originalTo = count.getToSequence();
                        if (needsGap) {
                            String from = count.getFromSequence();
                            counts.put("-", BaseInformationRecords.CountInfo.newBuilder().setFromSequence(from)
                                    .setToSequence("-").setMatchesReference(from.charAt(0) == '-')
                                    .setGenotypeCountForwardStrand(0).setGenotypeCountReverseStrand(0));
                            needsGap = false;
                        }
                        if (count.getToSequence().length() == 1) {

                            counts.put(count.getToSequence(), count);
                        } else {
                            if (count.getIsIndel()) {
                                // count is an indel count.
                                String adjustedTo = count.getToSequence();
                                if (adjustedTo.length() > offset) {
                                    adjustedTo = adjustedTo.substring(offset, offset + 1);
                                } else {
                                    LOG.warn(String.format("offset %d outside of to sequence %s.", offset, count.getToSequence()));
                                    adjustedTo = "-";
                                }
                                BaseInformationRecords.CountInfo.Builder countForBase = counts.get(adjustedTo);

                                // add the count the count of the indel to the count of the base:
                                countForBase.setGenotypeCountForwardStrand(countForBase.getGenotypeCountForwardStrand() + count.getGenotypeCountForwardStrand());
                                countForBase.setGenotypeCountReverseStrand(countForBase.getGenotypeCountReverseStrand() + count.getGenotypeCountReverseStrand());
                            }
                        }
                        countIndex++;
                    }
                    // save counts back in copy:
                    sample.clearCounts();
                    sample.addCounts(counts.get("A"));
                    sample.addCounts(counts.get("C"));
                    sample.addCounts(counts.get("T"));
                    sample.addCounts(counts.get("G"));
                    sample.addCounts(counts.get("-"));
                    sample.addCounts(counts.get("N"));
                    copy.setSamples(sampleIndex, sample);
                    sampleIndex++;
                }
                return copy;
            }


        }
    }
}
