package org.campagnelab.dl.genotype.segments;

import it.unimi.dsi.fastutil.ints.IntArrayList;
import it.unimi.dsi.fastutil.ints.IntList;
import it.unimi.dsi.fastutil.objects.*;
import org.campagnelab.dl.genotype.helpers.GenotypeHelper;
import org.campagnelab.dl.varanalysis.protobuf.BaseInformationRecords;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.util.*;
import java.util.function.Predicate;
import java.util.stream.Collectors;

public class RecordList implements Iterable<BaseInformationRecords.BaseInformation> {
    static private Logger LOG = LoggerFactory.getLogger(RecordList.class);

    @Override
    public Spliterator<BaseInformationRecords.BaseInformation> spliterator() {
        return records.spliterator();
    }

    ArrayList<BaseInformationRecords.BaseInformation> records = new ArrayList<>();
    /*
     * A map that contains list of records following a specific record. Use to combine indel position that are interleaved with
     * original records. When we need the full list of record, we combine record a in records with the records following a in
     * afterRecord.
     */
    public Object2ObjectOpenHashMap<BaseInformationRecords.BaseInformation, List<BaseInformationRecords.BaseInformation>> afterRecord = new Object2ObjectOpenHashMap<>();

    /**
     * Returns an iterator over elements of type {@code T}.
     *
     * @return an Iterator.
     */
    @Override
    public Iterator<BaseInformationRecords.BaseInformation> iterator() {
        return records.iterator();
    }

    /**
     * Get the first record, if there is one, or null otherwise.
     *
     * @return
     */
    public BaseInformationRecords.BaseInformation first() {
        if (records.isEmpty()) return null;
        return records.get(0);
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
        final BaseInformationRecords.BaseInformation builtCopy = copy.build();
        //  System.out.println("Adding builtCopy:"+FormatterCountHelper.format(builtCopy.getSamples(0)));
        addToFollowing(previous, builtCopy);
        //records.add(records.indexOf(previous) + 1, builtCopy);
        return previous;

    }

    void addToFollowing(BaseInformationRecords.BaseInformation previous, BaseInformationRecords.BaseInformation builtCopy) {
        List<BaseInformationRecords.BaseInformation> list = afterRecord.getOrDefault(previous, new ObjectArrayList<>());
        list.add(builtCopy);
        afterRecord.put(previous, list);
    }

    public void removeWhere(Predicate<BaseInformationRecords.BaseInformation> predicateIsTrue) {
        records.removeIf(predicateIsTrue);

    }

    public ObjectSet<BaseInformationRecords.BaseInformation> hideSet = new ObjectOpenHashSet<>();

    public void hideRecord(BaseInformationRecords.BaseInformation record) {
        hideSet.add(record);
    }
}
