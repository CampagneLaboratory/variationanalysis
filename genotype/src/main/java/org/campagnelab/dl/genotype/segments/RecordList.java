package org.campagnelab.dl.genotype.segments;

import it.unimi.dsi.fastutil.objects.Object2ObjectArrayMap;
import it.unimi.dsi.fastutil.objects.ObjectArrayList;
import org.campagnelab.dl.varanalysis.protobuf.BaseInformationRecords;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.util.Iterator;
import java.util.Map;

public class RecordList implements Iterable<BaseInformationRecords.BaseInformation> {
    static private Logger LOG = LoggerFactory.getLogger(RecordList.class);

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
