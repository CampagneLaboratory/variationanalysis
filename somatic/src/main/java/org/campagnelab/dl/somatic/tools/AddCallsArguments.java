package org.campagnelab.dl.somatic.tools;

import com.beust.jcommander.Parameter;
import com.beust.jcommander.Parameters;
import org.campagnelab.dl.framework.tools.arguments.ToolArguments;

/**
 * Arguments for AddCalls.
 */
@Parameters(commandDescription = "Add calls from mapped vcf to sbi/sbip files.")

public class AddCallsArguments implements ToolArguments {
    @Parameter(required = true, names = {"-i", "--input-files"}, description = "Input files in .bsi/.bsip format.")
    public String inputFile;

    @Parameter(required = true, names = {"-g", "--genotype-map"}, description = "Genotype may should have been generated with Goby's VCFToMapMode.")
    public String genotypeMap;

}

