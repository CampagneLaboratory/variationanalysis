package org.campagnelab.dl.genotype.tools;

import com.beust.jcommander.Parameter;
import com.beust.jcommander.Parameters;
import org.campagnelab.dl.framework.tools.arguments.ToolArguments;

/**
 * Arguments for AddTrueGenotypes.
 */
@Parameters(commandDescription = "Add calls from mapped vcf to sbi/sbip files.")

public class DebugGenotypeArguments implements ToolArguments {
    @Parameter(names = {"-i", "--input-filename"}, description = "Input files in .bsi/.bsip format.")
    public String inputFile;

    @Parameter(required = true, names = {"-m", "--genotype-map"}, description = "Genotype may should have been generated with Goby's VCFToMapMode.")
    public String genotypeMap;

    @Parameter(required = true, names = {"-g", "--genome"}, description = "Genome location to add calls with")
    public String genomeFilename;
    

}

