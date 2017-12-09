package org.campagnelab.dl.genotype.tools;

import com.beust.jcommander.Parameter;
import org.campagnelab.dl.framework.tools.arguments.ToolArguments;

/**
 * Arguments for the {@link FilterSSI} tool.
 *
 * @author manuele
 */
public class FilterSSIArguments implements ToolArguments {

    @Parameter(required = true, names = {"-i", "--input-file"}, variableArity = true,
            description = "Input SSI to filter, must be provided in .ssi format.")
    public String inputFile = null;

    @Parameter(required = true, names = {"-o", "--output-file"}, variableArity = true,
            description = "Output SSI with the filtered records.")
    public String outputFile = null;

    @Parameter(required = true, names = {"-s", "--start-position"}, variableArity = true,
            description = "Position from which the filter operates, must be provided in chr:pos format.")
    public String startPosition = null;


    @Parameter(required = true, names = {"-s", "--start-position"}, variableArity = true,
            description = "Position at which the filter stops, must be provided in chr:pos format.")
    public String endPosition = null;

}
