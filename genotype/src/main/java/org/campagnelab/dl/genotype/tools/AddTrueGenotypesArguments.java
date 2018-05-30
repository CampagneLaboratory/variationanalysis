package org.campagnelab.dl.genotype.tools;

import com.beust.jcommander.Parameter;
import com.beust.jcommander.Parameters;
import org.campagnelab.dl.framework.tools.arguments.ToolArguments;

/**
 * Arguments for AddTrueGenotypes.
 */
@Parameters(commandDescription = "Add calls from mapped vcf to sbi/sbip files.")

public class AddTrueGenotypesArguments implements ToolArguments {
    @Parameter(required = true, names = {"-i", "--input-filename"}, description = "Input files in .bsi/.bsip format.")
    public String inputFile;

    @Parameter(required = true, names = {"-o", "--output-filename"}, description = "Output filename for annoated data.")
    public String outputFilename;

    @Parameter(required = true, names = {"-m", "--genotype-map"}, description = "Genotype may should have been generated with Goby's VCFToMapMode.")
    public String genotypeMap;

    @Parameter(required = true, names = {"-g", "--genome"}, description = "Genome location to add calls with")
    public String genomeFilename;

    @Parameter(required = true, names = {"--confidence-regions"}, description = "A bed file with confidence regions. When provided, only sites within confidence regions are output.")
    public String confidenceRegionsFilename;

    @Parameter(required = false, names = {"-s", "--sample-index"}, description = "Add calls to an alternative sample in the sbi file (default if first sample, index 0")
    public int sampleIndex = 0;

    @Parameter( names = { "--ref-sampling-rate"}, description = "Sampling rate for positions where the true genotype matches the reference.")
    public float referenceSamplingRate=1.0f;

    @Parameter( names = { "--consider-indels"}, description = "When true, add true genotypes for indels. False (default) ignores indels.")
    public boolean considerIndels;

    @Parameter( names = { "--indels-as-ref"}, description = "When true, treat add the first base of indels as ref if they aren't considered. Ignored if indels considered.")
    public boolean indelsAsRef = true;
}

