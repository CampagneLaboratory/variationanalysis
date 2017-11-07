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
    @Parameter(  names = {"-l", "--level"}, description = "Simplify up to level (0,1,2).")
    public int simplifyLevel=3;

    @Parameter(  names = {"-I", "--indels"}, description = "Print only indels.")
    public boolean indels;
}

