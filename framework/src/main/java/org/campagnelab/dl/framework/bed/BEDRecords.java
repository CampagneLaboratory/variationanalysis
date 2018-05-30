package org.campagnelab.dl.framework.bed;

import it.unimi.dsi.fastutil.objects.Object2ObjectAVLTreeMap;
import it.unimi.dsi.fastutil.objects.Object2ObjectMap;
import it.unimi.dsi.fastutil.objects.ObjectArrayList;

import java.util.*;

/**
 * A collection of BED records.
 */
public class BEDRecords {
    final Object2ObjectMap<String, ObjectArrayList<BEDRecord>> store = new Object2ObjectAVLTreeMap<String, ObjectArrayList<BEDRecord>>();

    /**
     * Add a record to the set of records.
     *
     * @param record
     */
    public void add(BEDRecord record) {
        ObjectArrayList<BEDRecord> recordList = store.get(record.chromosome);
        if (recordList == null) {
            recordList = new ObjectArrayList<>();
            store.put(record.chromosome, recordList);
        }

        recordList.add(record);
    }

    /**
     * Sort records inside each chromosome by start position.
     */
    public void sort() {
        for (String chromosome: store.keySet()){
            List<BEDRecord> recordList = store.get(chromosome);
            recordList.sort(Comparator.comparingInt(value -> value.startPosition));
        }
    }
    Comparator<BEDRecord> COMPARATOR = (r1, r2) -> Integer.compare(r1.startPosition,r2.startPosition);
    /**
     * Returns true when the records overlap the specified range.
     * @param chromosome of the range.
     * @param startOfRange start of the range.
     * @param endOfRange end of the range.
     * @return
     */
    public boolean overlaps(String chromosome, int startOfRange, int endOfRange) {
        ObjectArrayList<BEDRecord> recordList = store.get(chromosome);
        if (recordList==null) {
            // unable to find chromosomes in regions
            return false;
        }
        BEDRecord key = new BEDRecord(chromosome, startOfRange, endOfRange);
        int r=Collections.binarySearch(recordList,key,COMPARATOR);
        if (r>=0) {
            // at least one record has the same start, so there is overlap.
            return true;
        }else {
            // calculate the insertion point: where the key would be if it was in the list.
            int i=-1-r;

            // now we need to look at surrounding elements to check for overlap:
            if (i-1>0 && recordList.get(i-1).overlapRange(startOfRange,endOfRange)) {
                return true;
            }
            if (i<recordList.size() &&recordList.get(i).overlapRange(startOfRange,endOfRange) ) {
                return true;
            }
            return false;
        }
    }

    public int numRecords() {
        Optional<Integer> result = this.store.values().stream().map(bedRecords -> bedRecords.size()).reduce((size1, size2) -> size1 + size2);
        return result.get();
    }

}

