package org.campagnelab.dl.genotype.segments.splitting;

import org.campagnelab.dl.genotype.segments.Segment;

import java.util.List;

/**
 * A strategy to split a segment.
 *
 * @author manuele
 */
public interface SplitStrategy{

    /**
     * Applies this strategy to the given argument.
     *
     * @param segment the function argument
     * @return the segment resulting from the splitting.
     */
    List<Segment> apply(Segment segment);
}
