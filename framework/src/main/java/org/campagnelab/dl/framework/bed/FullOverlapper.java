package org.campagnelab.dl.framework.bed;

/**
 * This class reports that every position overlaps with the record.
 */
public class FullOverlapper extends BEDRecords {
    @Override
    public boolean overlaps(String chromosome, int startOfRange, int endOfRange) {
        return true;
    }
}
