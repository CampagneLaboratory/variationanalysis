package org.campagnelab.dl.genotype.segments;

import java.util.function.Function;

/**
 * A function used to post-process a segment to rearrange indel genotypes as needed
 * to make a sequence with one base time step.
 */
public interface PostProcessSegmentFunction extends Function<Segment, Segment> {
}
