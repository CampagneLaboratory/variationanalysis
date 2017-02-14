package org.campagnelab.dl.genotype.tools;

import com.beust.jcommander.Parameter;
import com.beust.jcommander.Parameters;
import org.campagnelab.dl.framework.tools.arguments.ToolArguments;

/**
 * Arguments for SbiStats.
 */
@Parameters(commandDescription = "Get call stats from sbi/sbip files.")

public class SbiStatsArguments implements ToolArguments {
    @Parameter(required = true, names = {"-i", "--input-filename"}, description = "Input files in .bsi/.bsip format.")
    public String inputFile;


    @Parameter(required = false, names = {"-s", "--sample-index"}, description = "Add calls to an alternative sample in the sbi file (default if first sample, index 0")
    public int sampleIndex = 0;


}

