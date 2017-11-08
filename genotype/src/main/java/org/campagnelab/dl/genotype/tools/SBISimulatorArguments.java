package org.campagnelab.dl.genotype.tools;

import com.beust.jcommander.Parameter;
import org.campagnelab.dl.framework.tools.arguments.ToolArguments;

/**
 * Created by mas2182 on 11/3/17.
 */
public class SBISimulatorArguments implements ToolArguments {

    @Parameter(required = true, names = {"-o", "--output-filename"}, description = "Output filename for the sbi.")
    public String outputFilename;

    @Parameter(required = true, names = {"-i", "--input-file"}, variableArity = true,
            description = "Input file with the variant map, must be provided in .varmap format (produced with Goby3).")
    public String inputFile = null;

    @Parameter(required = false, names = {"-c", "--chromosome"}, variableArity = true,
            description = "Write sbi only for variants on the specified chromosome.")
    public String chromosome = null;

    @Parameter(names = "--read-N", description = "Read at most N chromosomes from the varmap, then stop.")
    public long readN=Long.MAX_VALUE;

    @Parameter(names ={"-v", "--verbose"}, description = "Be more verbose.")
    boolean verbose;
}

