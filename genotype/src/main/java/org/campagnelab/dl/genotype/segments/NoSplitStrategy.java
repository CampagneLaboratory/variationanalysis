package org.campagnelab.dl.genotype.segments;

import java.util.Arrays;
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
        return Arrays.asList(segment);
    }
}
