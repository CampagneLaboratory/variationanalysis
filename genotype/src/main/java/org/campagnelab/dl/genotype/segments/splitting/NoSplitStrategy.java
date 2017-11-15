package org.campagnelab.dl.genotype.segments.splitting;

import org.campagnelab.dl.genotype.segments.Segment;

import java.util.Collections;
import java.util.List;

/**
 * A no strategy. It does not split the segment.
 *
 */
public class NoSplitStrategy implements SplitStrategy {


    /**
     * Applies this strategy to the given argument.
     *
     * @param segment the function argument
     * @return the segment resulting from the splitting.
     */
    @Override
    public List<Segment> apply(Segment segment) {
        return Collections.singletonList(segment);
    }
}
