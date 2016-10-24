package org.campagnelab.dl.varanalysis.tools;

import com.beust.jcommander.Parameter;
import com.beust.jcommander.Parameters;

/**
 * Arguments for the Show tool.
 */

@Parameters(commandDescription = "Show records of a dataset. The dataset must be provided in the .sbi/sbip format.")

public class ShowArguments implements ToolArguments {
    @Parameter(required = true, names = {"-i", "--dataset"}, description = "Path to the dataset to show.")
    public String datasetFilename;

    @Parameter(names = {"-n", "--num-examples"}, description = "number of examples to show.")
    public int showN = Integer.MAX_VALUE;

}
