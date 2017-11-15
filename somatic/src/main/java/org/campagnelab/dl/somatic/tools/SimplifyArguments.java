package org.campagnelab.dl.somatic.tools;

import com.beust.jcommander.Parameter;
import com.beust.jcommander.Parameters;
import org.campagnelab.dl.framework.tools.arguments.ToolArguments;

/**
 * Arguments of Simplify tool.
 */
@Parameters(commandDescription = "Simplify an sbi/sbip file by dropping fields.")

public class SimplifyArguments implements ToolArguments {
    @Parameter(required = true,  names = {"-i", "--input-file"}, description = "Input file in .bsi/.bsip format.")
    public String inputFile;


    @Parameter(required = true, names = {"-o", "--output"}, description = "Output basename.")
    public String outputFile;

    @Parameter(required=false, names = {"-c", "--complexity"}, description = "Level of complexity left in the output. " +
            "Zero removes most data, keeping only very short messages, 1 produces a more complex output by removing less" +
            "fields, and so on, until 10 when the output is of similar complexity to the output, but with some fields" +
            "cleared.")
    int outputComplexityLevel=0;
}

