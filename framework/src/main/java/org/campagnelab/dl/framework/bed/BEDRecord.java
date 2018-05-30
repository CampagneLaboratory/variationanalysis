package org.campagnelab.dl.framework.bed;

/**
 * Represents a BED record, i.e., a segment/range on a chromosome.
 */
public class BEDRecord {
    /**
     * zero-based start position. The start position is included in the range, at the very start of it.
     */
    int startPosition;
    /**
     * zero-based end position (the end position is not part of the record, but just to the right of it).
     */
    int endPosition;
    /**
     * Chromosome where the segment/record is located.
     */
    public String chromosome;

    /**
     * Construct a new record.
     *
     * @param chromosome
     * @param segmentStart
     * @param segmentEnd
     */
    public BEDRecord(String chromosome, int segmentStart, int segmentEnd) {
        this.chromosome = chromosome;
        this.startPosition = segmentStart;
        this.endPosition = segmentEnd;

    }


    /**
     * Determine if the record overlaps with a single base. This method evaluates to position>=start && position<=end;
     *
     * @param position query position.
     * @return True or false.
     */
    public final boolean overlapPosition(final int position) {
        return position >= startPosition && position <= endPosition;
    }

    /**
     * Determine if the record overlaps a range.
     * @param startRange start of the range to test overlap with.
     * @param endRange end of the range to test overlap with.
     * @return true if the record overlaps the range, false otherwise.
     */
    public boolean overlapRange(final int startRange, final int endRange) {
        boolean overlap = false;

        //4-cases to consider

        //1st: segment2's start is between this segment start and end AND segment2's end is after this end
        if ((this.startPosition <= startRange && startRange <= this.startPosition) &&
                (endRange >= this.endPosition)) {
            overlap = true;
        }
        //2nd: this start is between segment2 start and end AND this end is after segment2 end
        else if ((startRange <= this.startPosition && this.startPosition <= endRange) &&
                (this.endPosition >= endRange)) {
            overlap = true;

        }
        //3rd: this is contained within segment2
        else if ((startRange <= this.startPosition && this.startPosition <= endRange) &&
                (startRange <= this.endPosition && this.endPosition <= endRange)) {
            overlap = true;

        }
        //4th: segment2 is contained within this
        else if ((this.startPosition <= startRange && startRange <= this.endPosition) &&
                (this.startPosition <= endRange && endRange <= this.endPosition)) {
            overlap = true;
        }

        return overlap;
    }
}
