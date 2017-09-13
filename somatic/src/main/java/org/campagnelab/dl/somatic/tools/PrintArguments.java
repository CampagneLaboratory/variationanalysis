package org.campagnelab.dl.somatic.tools;

import com.beust.jcommander.Parameter;
import com.beust.jcommander.Parameters;
import org.campagnelab.dl.framework.tools.arguments.ToolArguments;

/**
 * Arguments of Print tool.
 */
@Parameters(commandDescription = "Print an sbi file to protobuf text format.")

public class PrintArguments implements ToolArguments {
    @Parameter(required = true,  names = {"-i", "--input-file"}, description = "Input file in .bsi/.bsip format.")
    public String inputFile;

}

