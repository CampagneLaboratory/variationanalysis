package org.campagnelab.dl.genotype.tools;

import com.beust.jcommander.Parameter;
import com.beust.jcommander.Parameters;
import org.campagnelab.dl.framework.tools.arguments.ToolArguments;

/**
 * Arguments for AddTrueGenotypes.
 */
@Parameters(commandDescription = "Add calls from mapped vcf to sbi/sbip files.")

public class FilterSBIArguments implements ToolArguments {
    @Parameter(required = true, names = {"-i", "--input-filename"}, description = "Input files in .sbi/.sbip format.")
    public String inputFile;

    @Parameter(required = true, names = {"-o", "--output-filename"}, description = "Output filename for annoated data.")
    public String outputFilename;

    @Parameter(required = false, names = {"-s", "--sample-index"}, description = "Remove snps to an alternative sample in the sbi file (default if first sample, index 0")
    public int sampleIndex = 0;

    @Parameter( names = { "--other-sampling-rate"}, description = "Keep sites that are neither variants nor indels.true genotype matches the reference.")
    public double otherSamplingRate=0.001;

    @Parameter( names = { "seed"}, description = "optional custom random seed.")
    public int seed=240965;
    @Parameter(names={"--remove-SNPs"},description = "Remove SNPs when specified. Keep indels or reference matching sites.")
    public boolean removeSNPs;
    @Parameter(names={"--remove-ref-matching"},description = "Remove reference matching sites when specified. Keep SNPs or indels.")
    public boolean removeReferenceMatching;
}

