package org.campagnelab.dl.somatic.tools;

import com.beust.jcommander.Parameter;
import com.beust.jcommander.Parameters;
import org.campagnelab.dl.framework.tools.arguments.ToolArguments;
import org.campagnelab.dl.somatic.intermediaries.OneSampleCanonicalSimulationStrategy;
import org.campagnelab.dl.somatic.intermediaries.TwoSampleCanonicalSimulationStrategy;

import java.util.List;

/**
 * Created by rct66 on 11/14/16.
 * jcommander arguments for mutator2
 */
@Parameters(commandDescription = "Arguments for planting mutations in an sbi file.")

public class Mutator2Arguments implements ToolArguments {
    @Parameter(required = true, names = {"-i", "--input-file"}, description = "Input file in .bsi/.bsip format.")
    public String inputFile;

    @Parameter(required = false, names = {"-c", "--canonical-threshold"}, description = "Fraction of total counts used to determine whether both samples have matching, canonical zygosity.")
    public Float canonThreshold = 0.9f;

    @Parameter(required = false, names = {"-f", "--zygosity-floor"}, description = "Fraction of total counts the base with the second-most counts must have for the record to be considered heterozygous")
    public Float heteroHeuristic = 0.1f;

    @Parameter(required = true, names = {"-o", "--output-file"}, description = "Path to output file.")
    public String outputFile;

    @Parameter(required = true, names = {"-s", "--strategy"}, description = "Strategy classname.")
    public String strategyClassname= OneSampleCanonicalSimulationStrategy.class.getCanonicalName();

    @Parameter( names = { "--random-seed"}, description = "Random seed.")
    public long seed= 2398823;
}
