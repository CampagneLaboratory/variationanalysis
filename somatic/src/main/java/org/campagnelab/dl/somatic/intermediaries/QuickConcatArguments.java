package org.campagnelab.dl.somatic.intermediaries;

import com.beust.jcommander.Parameter;
import org.campagnelab.dl.framework.tools.arguments.ToolArguments;

import java.util.ArrayList;
import java.util.List;

/**
 * Created by fac2003 on 10/17/16.
 */
public class QuickConcatArguments implements ToolArguments {
    @Parameter(required = true, variableArity = true, names = {"-i", "--input-files"}, description = "Input files in .bsi/.bsip format.")
    public List<String> inputFiles = new ArrayList<>();

     @Parameter(required = true, names = {"-o", "--output-prefix"}, description = "Prefix for the output filenames.")
    public String outputFile;

    @Parameter(required = false, names = {"-b", "--buffer-size"}, description = "Size of the copy buffer. Should be large enough to be efficient, but small enough that quick concat can be interrupted by ^C. ")
    public int copyBufferSize = 100*1024*1024;

    @Parameter(names={"-f","--force"}, description = "Force override the output if it already exists.")
    public    boolean force=false;
}
