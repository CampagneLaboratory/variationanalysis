package org.campagnelab.dl.somatic.tools;

import com.beust.jcommander.Parameter;
import com.beust.jcommander.Parameters;
import org.campagnelab.dl.framework.tools.arguments.ToolArguments;

import java.util.List;

/**
 * Created by rct2002 on 10/23/17.
 */
@Parameters(commandDescription = "Split a file into several components by matching chromosome name. -v goes to validation, -t goes to test, and non-matches go to training.")

public class SplitByChromosomeArguments implements ToolArguments {
    @Parameter(required = true, names = {"-i", "--input-file"}, description = "Input file in .bsi/.bsip format.")
    public String inputFile;

    @Parameter(required = true, names = {"-v", "--validation"}, description = "set of chromosomes to include in the validation set, " +
            "for example -v chr20 -v chr21")
    List<String> valChromosomes;

    @Parameter(required = true, names = {"-t", "--test"}, description = "set of chromsomes to include in the test set " +
            "for example -v chr22 -v chr23")
    List<String> testChromosomes;

    @Parameter(required=true, names = {"-o", "--output-prefix"}, description = "Prefix for the output filenames.")
    public String outputFile;

}

