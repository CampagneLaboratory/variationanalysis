package org.campagnelab.dl.genotype.segments;

import it.unimi.dsi.fastutil.ints.Int2ObjectMap;
import it.unimi.dsi.fastutil.ints.Int2ObjectOpenHashMap;
import it.unimi.dsi.fastutil.objects.*;
import it.unimi.dsi.lang.MutableString;
import org.campagnelab.dl.genotype.helpers.GenotypeHelper;
import org.campagnelab.dl.varanalysis.protobuf.BaseInformationRecords;
import org.campagnelab.dl.varanalysis.protobuf.SegmentInformationRecords;
import org.campagnelab.goby.baseinfo.SequenceSegmentInformationWriter;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Spliterator;
import java.util.function.Consumer;
import java.util.function.Function;
import java.util.function.Predicate;
import java.util.stream.Collectors;
import java.util.stream.Stream;

/**
 * Holds the current open segment before it is stored in the list.
 */
public class Segment {
    public final Function<BaseInformationRecords.BaseInformation, SegmentInformationRecords.Base.Builder> fillInFeatures;
    private int firstPosition = 0;
    private int firstReferenceIndex = 0;
    private String firstReferenceId = "";
    private int lastPosition = 0;
    private String lastReferenceId = "";
    private int lastReferenceIndex = 0;
    public RecordList recordList = new RecordList();

    public Segment(Function<BaseInformationRecords.BaseInformation, SegmentInformationRecords.Base.Builder> fillInFeatures, BaseInformationRecords.BaseInformation first) {
        //System.out.println("Open a new segment at ref " + first.getReferenceId() + " position " + Integer.toString(first.getPosition()));
        this.add(first);
        this.firstPosition = first.getPosition();
        this.firstReferenceIndex = first.getReferenceIndex();
        this.firstReferenceId = first.getReferenceId();
        this.fillInFeatures = fillInFeatures;
    }

    private void setAsLast(BaseInformationRecords.BaseInformation record) {
        this.lastPosition = record.getPosition();
        this.lastReferenceId = record.getReferenceId();
        this.lastReferenceIndex = record.getReferenceIndex();
    }

    /**
     * Write the segment to an SBI writer.
     *
     * @param segmentConsumer
     */
    public void construct(Consumer<SegmentInformationRecords.SegmentInformation> segmentConsumer) {
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
        int numSamples = recordList.first().getSamplesCount();
        SegmentInformationRecords.Sample.Builder sampleBuilder[] = new SegmentInformationRecords.Sample.Builder[numSamples];
        for (int sampleIndex = 0; sampleIndex < numSamples; sampleIndex++) {

            sampleBuilder[sampleIndex] = SegmentInformationRecords.Sample.newBuilder();

        }

        getAllRecords().forEach(record -> {

            for (int sampleIndex = 0; sampleIndex < numSamples; sampleIndex++) {
                SegmentInformationRecords.Base.Builder base = fillInFeatures.apply(record);

                sampleBuilder[sampleIndex].addBase(base);

                //System.out.println("New base " + segmentStats[0] );


            }
        });
        for (int sampleIndex = 0; sampleIndex < numSamples; sampleIndex++) {
            builder.addSample(sampleBuilder[sampleIndex]);
        }

        final SegmentInformationRecords.SegmentInformation built = builder.build();

        segmentConsumer.accept(built);

    }

    /**
     * Adds a record to the current segment
     *
     * @param record
     */
    public void add(BaseInformationRecords.BaseInformation record) {
        this.recordList.add(record);
        this.setAsLast(record);
    }

    @Override
    public String toString() {
        return String.format("Segment{start=%s:%d end=%s:%d length=%d}%n", getFirstReferenceId(), getFirstPosition(),
                getLastReferenceId(), getLastPosition(), actualLength());
    }

    public int actualLength() {
        if (recordList.size() == 0)
            return 0;
        else
            return recordList.size() + recordList.afterRecord.size();
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

    private Int2ObjectMap positionsToTrueGenotypes = new Int2ObjectOpenHashMap();

    public void populateTrueGenotypes() {
        positionsToTrueGenotypes.clear();
        for (BaseInformationRecords.BaseInformation record : recordList) {
            String trueGenotype = record.getTrueGenotype();
            if (trueGenotype.length() == 3) {
                // no indel, simple A/B genotype:
                positionsToTrueGenotypes.put(record.getPosition(), trueGenotype);
            }
        }
    }

    public String getTrueGenotype(int position) {
        return (String) positionsToTrueGenotypes.getOrDefault(position, "-");
    }

    /**
     * Insert a copy after a record. The copy has the same position as the record, but follows in order (usually
     * representing a position contributing to an indel).
     *
     * @param record after which copy will be inserted
     * @param copy   to insert.
     */
    public void insertAfter(BaseInformationRecords.BaseInformation record, BaseInformationRecords.BaseInformation.Builder copy) {

        recordList.addToFollowing(record, copy.build());
    }

    /**
     * Returns the complete list of records, including those interleaved with genomic positions (for insertion/deletion).
     *
     * @return
     */
    public Iterable<BaseInformationRecords.BaseInformation> getAllRecords() {
        ObjectArrayList<BaseInformationRecords.BaseInformation> list = new ObjectArrayList(recordList.size() * 3 / 2);
        for (BaseInformationRecords.BaseInformation record : recordList) {
            if (!recordList.hideSet.contains(record)) {
                list.add(record);
            }
            list.addAll(recordList.afterRecord.getOrDefault(record, Collections.emptyList()));
        }

        return list;
    }

    public Spliterator<BaseInformationRecords.BaseInformation> getAllRecordsSplit() {
        ArrayList<BaseInformationRecords.BaseInformation> list = new ArrayList<>();
        for (BaseInformationRecords.BaseInformation record : recordList) {
            list.add(record);
            list.addAll(recordList.afterRecord.getOrDefault(record, Collections.emptyList()));
        }
        return list.parallelStream().spliterator();
    }

    public static String showGenotypes(SegmentInformationRecords.SegmentInformation segmentInformation) {
        MutableString result = new MutableString();
        // for (int sampleIndex = 0; sampleIndex < segmentInformation.getSampleCount(); sampleIndex++) {
        int sampleIndex = 0;
        for (SegmentInformationRecords.Base base : segmentInformation.getSample(sampleIndex).getBaseList()) {
            ArrayList<String> alleles = new ArrayList<>();
            for (String allele : base.getTrueLabelList()) {
                alleles.add(allele);
            }

            result.append(String.format("ref=%s\ttrueGenotype=%s\tcounts=%s%n", base.getReferenceAllele(),
                    alleles.stream().collect(Collectors.joining("/")),
                    base.getFormattedCounts()));
        }
        return result.toString();

    }

    public void remove(BaseInformationRecords.BaseInformation record) {
        recordList.records.remove(record);
    }

    public void removeWhere(Predicate<BaseInformationRecords.BaseInformation> predicateIsTrue) {
        recordList.records.removeIf(predicateIsTrue);
    }
}