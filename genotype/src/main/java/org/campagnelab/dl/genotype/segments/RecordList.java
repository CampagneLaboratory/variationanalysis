package org.campagnelab.dl.genotype.segments;

import it.unimi.dsi.fastutil.ints.IntArrayList;
import it.unimi.dsi.fastutil.objects.Object2ObjectOpenHashMap;
import it.unimi.dsi.fastutil.objects.ObjectArrayList;
import it.unimi.dsi.fastutil.objects.ObjectOpenHashSet;
import it.unimi.dsi.fastutil.objects.ObjectSet;
import org.campagnelab.dl.genotype.helpers.GenotypeHelper;
import org.campagnelab.dl.varanalysis.protobuf.BaseInformationRecords;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import java.util.Spliterator;
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
    Object2ObjectOpenHashMap<BaseInformationRecords.BaseInformation, List<BaseInformationRecords.BaseInformation>> afterRecord = new Object2ObjectOpenHashMap<>();

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

    /**
     * Adjust counts to reduce indels to a single base, using offset and the current indel sequences to determine where
     * to increase the base count.
     *
     * @param copy   record that needs to have counts adjusted.
     * @param offset offset inside the indel sequence, which identifies the base to increment.
     * @return a builder where the adjustment has been made.
     */
    public BaseInformationRecords.BaseInformation.Builder adjustCounts(BaseInformationRecords.BaseInformation.Builder copy, int offset) {
        // we store counts in a map for easy access (map keyed on to sequence of the count):
        // we know we may need a gap count, so we add one, because none in the sbi:

        int sampleIndex = 0;
        int countIndex = 0;
        String recordTrueGenotype = copy.getTrueGenotype();

        for (BaseInformationRecords.SampleInfo.Builder sample : copy.getSamplesBuilderList()) {
            ObjectArrayList<BaseInformationRecords.CountInfo.Builder> countsToKeep = new ObjectArrayList();

            for (BaseInformationRecords.CountInfo.Builder count : sample.getCountsBuilderList()) {

                if (count.getToSequence().length() != 1) {
                    if (count.getIsIndel()) {
                        // count is an indel count.
                        String adjustedTo = count.getToSequence();
                        String adjustedFrom = count.getFromSequence();
                        String adjustedTrueGenotype = recordTrueGenotype;
                        if (adjustedTo.length() > offset && adjustedFrom.length() > offset) {
                            adjustedFrom = adjustedFrom.substring(offset, offset + 1);
                            adjustedTo = adjustedTo.substring(offset, offset + 1);
                            adjustedTrueGenotype = adjusteTrueGenotype(adjustedTrueGenotype,offset);
                            count.setFromSequence(adjustedFrom);
                            copy.setReferenceBase(adjustedFrom);
                            count.setToSequence(adjustedTo);
                            countsToKeep.add(count);
                            sample.setTrueGenotype(adjustedTrueGenotype);
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
            for (BaseInformationRecords.CountInfo.Builder count : countsToKeep) {
                sample.addCounts(count);

            }
            copy.setSamples(sampleIndex, sample);
            sampleIndex++;

        }
        return copy;
    }

    private String adjusteTrueGenotype(String adjustedTrueGenotype, int offset) {

        ArrayList<String> alleles=new ArrayList();
        for (String allele : GenotypeHelper.getAlleles(adjustedTrueGenotype)) {
          if (allele.length()>offset) {
              alleles.add(allele.substring(offset, offset + 1));
          } else {
                alleles.add("-");
            }
        }
       return alleles.stream().collect(Collectors.joining("/"));

    }

    private void makeEmptyCount(ObjectArrayList<BaseInformationRecords.CountInfo.Builder> countsToKeep, BaseInformationRecords.CountInfo.Builder count) {
        count.setGenotypeCountForwardStrand(0);
        count.setGenotypeCountReverseStrand(0);
        countsToKeep.add(count);
    }


    public void removeWhere(Predicate<BaseInformationRecords.BaseInformation> predicateIsTrue) {
        records.removeIf(predicateIsTrue);

    }
    ObjectSet<BaseInformationRecords.BaseInformation> hideSet=new ObjectOpenHashSet<>();
    public void hideRecord(BaseInformationRecords.BaseInformation record) {
        hideSet.add(record);
    }
}
