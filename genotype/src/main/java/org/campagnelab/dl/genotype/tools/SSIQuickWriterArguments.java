package org.campagnelab.dl.genotype.tools;

import com.beust.jcommander.Parameter;
import com.beust.jcommander.Parameters;
import org.campagnelab.dl.framework.tools.arguments.ToolArguments;

/**
 *
 * Arguments for the {@link SSIQuickWriter} tool.
 * @author joshuacohen
 */
@Parameters(commandDescription = "Convert a set of Single Base Information data to Sequence Segment Information structures.")
public class SSIQuickWriterArguments implements ToolArguments {

    @Parameter(required = true, names = {"-o", "--outputFile"}, description = "Path to output SSI file to write to")
    public String outputFile;

    @Parameter(names = {"-n", "--num-segments"}, description = "Number of Segments to write")
    public int numSegments = 90;

    @Parameter(names = {"-b", "--bases-per-segment"}, description = "Number of Bases to write for each segment")
    public int basesPerSegment = 900;

    @Parameter(names = {"-f", "--features-per-base"}, description = "Number of features to write for each Base")
    public int featuresPerBase = 1800;

    @Parameter(names = {"-u", "--upper-bound"}, description = "Upper bound on value of features to generate")
    public int upperBound = 18;
}
