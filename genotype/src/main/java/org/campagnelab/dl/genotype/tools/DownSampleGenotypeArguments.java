package org.campagnelab.dl.genotype.tools;

import com.beust.jcommander.Parameter;
import com.beust.jcommander.Parameters;
import org.campagnelab.dl.framework.tools.arguments.ToolArguments;

/**
 * Arguments for DownSampleGenotypes.
 */
@Parameters(commandDescription = "Down-sample genotypes in an sbi/sbip file according to their type. Useful to produce " +
        "files more focused on indels/and or hets.")

public class DownSampleGenotypeArguments implements ToolArguments {
    @Parameter(required = true, names = {"-i", "--input-filename"}, description = "Input file in .bsi/.bsip format.")
    public String inputFile;

    @Parameter(required = true, names = {"-o", "--output-filename"}, description = "Output filename where down-sampled data will be written.")
    public String outputFilename;

    @Parameter(required = true, names = {"-g", "--genome"}, description = "Genome location to add calls with")
    public String genomeFilename;

    @Parameter(required = false, names = {"-s", "--sample-index"}, description = "Indicate the sample that should be down-sampled (default if first sample, index 0")
    public int sampleIndex = 0;

    @Parameter(required = false, names = { "--keep-indels"}, description = "Do not down-sample indel sites.")
    public boolean keepAllIndels;

    @Parameter(required = false, names = { "--keep-hets"}, description = "Do not down-sample heterozygous sites.")
    public boolean keepAllHeterozygotes;

    @Parameter( names = { "--other-sampling-rate"}, description = "Sampling rate for genotypes that match the down-sampling criteria. Default is not to include them.")
    public float otherSamplingRate =0f;

    @Parameter( names = { "--balancing-ratio"}, description = "Use 2 to downsample such that the number of non selected genotypes if about double the number of selected genotypes. " +
            "This argument overrides --other-sampling-rate by continuously adjusting the rate according to the number of selected genotypes found so far. Currently not implemented.")
    public Float balancingRatio =null;

    @Parameter( names = { "seed"}, description = "optional custom random seed.")
    public int seed=240965;

}

