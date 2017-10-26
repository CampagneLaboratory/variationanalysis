package org.campagnelab.dl.genotype.segments;

public class SnpOnlyPostProcessSegmentFunction implements PostProcessSegmentFunction {
    @Override
    public Segment apply(Segment segment) {

        // remove all indels
        segment.recordList.removeWhere(record -> record.getTrueGenotype().length() > 3);
        return segment;

    }
}
