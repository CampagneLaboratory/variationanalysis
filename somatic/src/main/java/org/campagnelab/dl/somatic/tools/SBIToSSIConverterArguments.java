package org.campagnelab.dl.somatic.tools;

import com.beust.jcommander.Parameter;
import com.beust.jcommander.Parameters;
import org.campagnelab.dl.framework.tools.arguments.ToolArguments;

import java.util.ArrayList;
import java.util.List;

/**
 *
 * Arguments for the {@link SBIToSSIConverter} tool.
 * @author manuele
 */
@Parameters(commandDescription = "Convert a set of Single Base Information data to Sequence Segment Information structures.")
public class SBIToSSIConverterArguments implements ToolArguments {

    @Parameter(required = true, names = {"-i", "--input-file"}, variableArity = true,
            description = "Input file with the SBI, must be provided in .sbi/.sbip format (produced with Goby3).")
    public String inputFile = null;

    @Parameter(names = {"-g", "--gap"}, description = "Gap between two segments The default is 1.")
    public int gap = 1;

    @Parameter(names = "--output-prefix", description = "Prefix for the output saved files. ")
    public String ssiPrefix = null;
}
