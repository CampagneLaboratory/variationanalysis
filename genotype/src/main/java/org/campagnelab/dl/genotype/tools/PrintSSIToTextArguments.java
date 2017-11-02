package org.campagnelab.dl.genotype.tools;

import com.beust.jcommander.Parameter;
import org.campagnelab.dl.framework.tools.arguments.ToolArguments;

/**
 * Arguments for the {@link PrintSSIToText} tool.
 *
 * @author manuele
 */
public class PrintSSIToTextArguments implements ToolArguments {

    @Parameter(required = true,  names = {"-i", "--input-file"}, description = "Input file in .ssi/.ssip format.")
    public String inputFile;

    @Parameter(names = { "--no-features"}, description = "Do not show features.")
    public boolean removeFeatures;

    @Parameter(names = { "--no-labels"}, description = "Do not show labels.")
    public boolean removeLabels;
}
