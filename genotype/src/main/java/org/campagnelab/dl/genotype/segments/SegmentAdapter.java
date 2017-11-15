package org.campagnelab.dl.genotype.segments;

import org.campagnelab.dl.somatic.storage.RecordReader;
import org.campagnelab.dl.varanalysis.protobuf.BaseInformationRecords;
import org.jetbrains.annotations.NotNull;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.util.Iterator;
import java.util.Spliterator;

/**
 * An adapter to create segments from SBI records. Uses a minimum gap between records to determine how to group
 * records into contiguous segments.
 */
public class SegmentAdapter implements Iterable<Segment> {
    static private Logger LOG = LoggerFactory.getLogger(SegmentAdapter.class);
    private final RecordReader reader;
    private final int gap;

    public SegmentAdapter(RecordReader reader, int gap) {
        this.reader = reader;
        this.gap=gap;
    }

    @NotNull
    @Override
    public Iterator<Segment> iterator() {
        return new SegmentIterator(this,gap);
    }


    @Override
    public Spliterator<Segment> spliterator() {
        return null;
    }

    private class SegmentIterator implements Iterator<Segment> {
        private SegmentAdapter adapter;
        private int gap = 50;
        private Segment currentSegment = null;
        private Segment result;

        public SegmentIterator(SegmentAdapter adapter, int gap) {
            this.adapter = adapter;
            this.gap=gap;
        }

        @Override
        public boolean hasNext() {

                BaseInformationRecords.BaseInformation record = reader.nextRecord();
                boolean newSegment = false;
                if (record != null) {

                    while (!newSegment) {
                        final boolean sameSegment = (record.getPosition() - currentSegment.getLastPosition() <= gap) &&
                                record.getReferenceIndex() == currentSegment.getLastReferenceIndex();
                        if (sameSegment) {
                            currentSegment.add(record);
                        } else {
                            newSegment = true;
                        }
                    }
                    return currentSegment != null;
                } else {
                    return false;
                }

        }

        @Override
        public Segment next() {
            result = currentSegment;
            currentSegment = null;
            return result;
        }
    }
}
