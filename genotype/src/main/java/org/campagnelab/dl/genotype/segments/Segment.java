package org.campagnelab.dl.genotype.segments;

import com.google.common.collect.Iterators;
import it.unimi.dsi.fastutil.objects.*;
import it.unimi.dsi.lang.MutableString;
import org.campagnelab.dl.genotype.helpers.GenotypeHelper;
import org.campagnelab.dl.varanalysis.protobuf.BaseInformationRecords;
import org.campagnelab.dl.varanalysis.protobuf.SegmentInformationRecords;

import java.util.*;
import java.util.function.Consumer;
import java.util.function.Function;
import java.util.function.Predicate;
import java.util.stream.Collectors;

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
    protected RecordList recordList = new RecordList();

    public boolean isIndicesAdded() {
        return indicesAdded;
    }

    public void setIndicesAdded(boolean indicesAdded) {
        this.indicesAdded = indicesAdded;
    }

    private boolean indicesAdded;

    public Segment(Function<BaseInformationRecords.BaseInformation, SegmentInformationRecords.Base.Builder> fillInFeatures, BaseInformationRecords.BaseInformation first) {
        this.add(first);
        this.firstPosition = first.getPosition();
        this.firstReferenceIndex = first.getReferenceIndex();
        this.firstReferenceId = first.getReferenceId();
        this.fillInFeatures = fillInFeatures;
    }

    protected Segment(Function<BaseInformationRecords.BaseInformation, SegmentInformationRecords.Base.Builder> fillInFeatures) {
        this.fillInFeatures = fillInFeatures;
    }

    protected void setAsLast(BaseInformationRecords.BaseInformation record) {
        this.lastPosition = record.getPosition();
        this.lastReferenceId = record.getReferenceId();
        this.lastReferenceIndex = record.getReferenceIndex();
    }

    protected void setAsFirst(BaseInformationRecords.BaseInformation base) {
        this.firstPosition = base.getPosition();
        this.firstReferenceIndex = base.getReferenceIndex();
        this.firstReferenceId = base.getReferenceId();
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

        final long[] segmentStats = {0L, 0L, 0L};
        int numSamples = getFirstRecord().getSamplesCount();
        SegmentInformationRecords.Sample.Builder sampleBuilder[] = new SegmentInformationRecords.Sample.Builder[numSamples];
        for (int sampleIndex = 0; sampleIndex < numSamples; sampleIndex++) {
            sampleBuilder[sampleIndex] = SegmentInformationRecords.Sample.newBuilder();
        }
        getAllRecords().forEach(record -> {

            for (int sampleIndex = 0; sampleIndex < numSamples; sampleIndex++) {
                SegmentInformationRecords.Base.Builder base = fillInFeatures.apply(record);
                sampleBuilder[sampleIndex].addBase(base);
            }
        });
        builder.setLength(sampleBuilder[0].getBaseCount());
        for (int sampleIndex = 0; sampleIndex < numSamples; sampleIndex++) {
            builder.addSample(sampleBuilder[sampleIndex]);
        }
        final SegmentInformationRecords.SegmentInformation built = builder.build();
        segmentConsumer.accept(built);
    }

    /**
     * First record in the segment.
     *
     * @return
     */
    public BaseInformationRecords.BaseInformation getFirstRecord() {
        return recordList.first();
    }

    /**
     * Adds a record to the current segment
     *
     * @param record
     */
    public void add(BaseInformationRecords.BaseInformation record) {
        this.recordList.add(record);
        if (this.recordList.size() == 1)
            this.setAsFirst(record);
        this.setAsLast(record);
    }

    @Override
    public String toString() {
        return String.format("Segment{start=%s:%d end=%s:%d length=%d}%n", getFirstReferenceId(), getFirstPosition(),
                getLastReferenceId(), getLastPosition(), actualLength());
    }

    public int  actualLength() {
        if (recordList.size() == 0)
            return 0;
        else {
            //
            return Iterators.size(getAllRecords().iterator());
        }
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
    public Iterable<BaseInformationRecords.BaseInformation> getAllRecords(int startPosition, int endPosition) {
        ObjectArrayList<BaseInformationRecords.BaseInformation> list = new ObjectArrayList();
        for (BaseInformationRecords.BaseInformation record : recordList) {
            if (record.getPosition() >= startPosition && record.getPosition() < endPosition) {
                if (!getHiddenRecords().contains(record)) {
                    list.add(record);
                }
                list.addAll(getAfterRecords().getOrDefault(record, Collections.emptyList()));
            }
        }

        return list;
    }

    /**
     * Determine the length of this segment between start and end position, including the number of afterRecord
     * included in the range.
     *
     * @param startPosition start position, inclusive
     * @param endPosition   end position, inclusive
     * @return the number of records such that startPosition<=record.position=endPosition
     */
    public int actualLength(int startPosition, int endPosition) {
        int count = 0;
        for (BaseInformationRecords.BaseInformation record : recordList) {
            if (record.getPosition() >= startPosition && record.getPosition() <= endPosition) {
                if (!getHiddenRecords().contains(record)) {
                    count += 1;
                }
                count += getAfterRecords().getOrDefault(record, Collections.emptyList()).size();
            }
        }
        return count;
    }

    private ObjectArrayList<BaseInformationRecords.CountInfo.Builder> countsToKeep = new ObjectArrayList();

    /**
     * Adjust counts to reduce indels to a single base, using offset and the current indel sequences to determine where
     * to increase the base count.
     *
     * @param copy   record that needs to have counts adjusted.
     * @param offset offset inside the indel sequence, which identifies the base to increment.
     * @param longestReference
     * @return a builder where the adjustment has been made.
     */
    public BaseInformationRecords.BaseInformation.Builder adjustCounts(BaseInformationRecords.BaseInformation.Builder copy, int offset, String longestReference) {
        // we store counts in a map for easy access (map keyed on to sequence of the count):
        // we know we may need a gap count, so we add one, because none in the sbi:

        int sampleIndex = 0;
        int countIndex = 0;
        String recordTrueGenotype = expand(copy.getTrueGenotype(),longestReference);

        for (BaseInformationRecords.SampleInfo.Builder sample : copy.getSamplesBuilderList()) {
            countsToKeep.clear();
            sample.setPrePostProcessingGenotype(recordTrueGenotype);

            for (BaseInformationRecords.CountInfo.Builder count : sample.getCountsBuilderList()) {
                count.setOffset(0);
                if (count.getToSequence().length() != 1) {
                    if (count.getIsIndel() || count.getToSequence().length() > 1) {

                        // count is an indel count.
                        String adjustedTo = count.getToSequence();
                        String adjustedFrom = count.getFromSequence();
                        String adjustedTrueGenotype = recordTrueGenotype;
                        if (adjustedTo.length() > offset && adjustedFrom.length() > offset) {
                            adjustedFrom = adjustedFrom.substring(offset, offset + 1);
                            adjustedTo = adjustedTo.substring(offset, offset + 1);
                            adjustedTrueGenotype = adjustTrueGenotype(adjustedTrueGenotype, offset);
                            count.setFromSequence(adjustedFrom);
                            copy.setReferenceBase(adjustedFrom);
                            count.setToSequence(adjustedTo);
                            countsToKeep.add(count);
                            sample.setTrueGenotype(adjustedTrueGenotype);
                            sample.setPrePostProcessingGenotype(recordTrueGenotype);
                            count.setOffset(offset);
                            copy.setTrueGenotype(adjustedTrueGenotype);
                        } else {
                            // this indel does not contribute to the counts at this offset:
                            makeEmptyCount(countsToKeep, count);
                        }
                    } else {
                        // this base does not occur at this offset:
                        makeEmptyCount(countsToKeep, count);
                    }

                } else {
                    makeEmptyCount(countsToKeep, count);
                }
            }
            sample.clearCounts();
            countsToKeep = aggregateCounts(countsToKeep);
            for (BaseInformationRecords.CountInfo.Builder count : countsToKeep) {
                sample.addCounts(count);

            }
            copy.setSamples(sampleIndex, sample);
            sampleIndex++;

        }
        return copy;
    }

    private String expand(String trueGenotype, String longestReference) {
        String result="";
        for (String allele: GenotypeHelper.getAlleles(trueGenotype)) {
            if (longestReference.startsWith(allele)) {
                allele=longestReference;
            }
            result+=allele+"/";
        }
        return result.substring(0, result.length()-1);
    }

    private String adjustTrueGenotype(String adjustedTrueGenotype, int offset) {

        ArrayList<String> alleles = new ArrayList();
        for (String allele : GenotypeHelper.getAlleles(adjustedTrueGenotype)) {
            if (allele.length() > offset) {
                alleles.add(allele.substring(offset, offset + 1));
            } else {
                alleles.add("N");
            }
        }
        return alleles.stream().collect(Collectors.joining("/"));

    }

    private void makeEmptyCount(ObjectArrayList<BaseInformationRecords.CountInfo.Builder> countsToKeep, BaseInformationRecords.CountInfo.Builder count) {
        count.setGenotypeCountForwardStrand(0);
        count.setGenotypeCountReverseStrand(0);
        countsToKeep.add(count);
    }

    Comparator<BaseInformationRecords.CountInfo.Builder> comparator = (o1, o2) -> {
        int diff = o1.getFromSequence().compareTo(o2.getFromSequence());
        if (diff != 0) return diff;
        diff = o1.getToSequence().compareTo(o2.getToSequence());
        if (diff != 0) return diff;
        return 0;
    };


    private ObjectSet<BaseInformationRecords.CountInfo.Builder> toRemove = new ObjectArraySet();

    // this method aggregates forward and reverse counts when the from/to matches across count instances.
    // it keeps only one instance, with the sum of counts over redundant instances.
    private ObjectArrayList<BaseInformationRecords.CountInfo.Builder> aggregateCounts(ObjectArrayList<BaseInformationRecords.CountInfo.Builder> countsToKeep) {

        toRemove.clear();
        Collections.sort(countsToKeep, comparator);


        BaseInformationRecords.CountInfo.Builder previous = null;
        for (BaseInformationRecords.CountInfo.Builder count : countsToKeep) {
            if (previous != null) {
                if (comparator.compare(previous, count) == 0) {
                    int sumForward = previous.getGenotypeCountForwardStrand() + count.getGenotypeCountForwardStrand();
                    int sumReverse = previous.getGenotypeCountReverseStrand() + count.getGenotypeCountReverseStrand();
                    previous.setGenotypeCountForwardStrand(sumForward);
                    previous.setGenotypeCountReverseStrand(sumReverse);
                    previous.setMatchesReference(previous.getFromSequence().equals(previous.getToSequence()));
                    toRemove.add(count);
                }
            }
            previous = count;

        }
        countsToKeep.removeAll(toRemove);
        return countsToKeep;
    }

    public Iterable<BaseInformationRecords.BaseInformation> getAllRecords() {
        return getAllRecords(-1, Integer.MAX_VALUE);
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

    public BaseInformationRecords.BaseInformation getRecordAt(int position) {
        for (BaseInformationRecords.BaseInformation record : getAllRecords()) {
            if (record.getPosition() == position)
                return record;
        }
        return null;
    }

    public void hideRecord(BaseInformationRecords.BaseInformation record) {
        this.recordList.hideRecord(record);
    }

    public ObjectSet<BaseInformationRecords.BaseInformation> getHiddenRecords() {
        return this.recordList.hideSet;
    }

    public Object2ObjectOpenHashMap<BaseInformationRecords.BaseInformation, List<BaseInformationRecords.BaseInformation>> getAfterRecords() {
        return this.recordList.afterRecord;
    }
}