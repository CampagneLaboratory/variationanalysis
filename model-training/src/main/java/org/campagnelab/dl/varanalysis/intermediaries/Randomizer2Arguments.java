package org.campagnelab.dl.varanalysis.intermediaries;

import com.beust.jcommander.Parameter;
import com.beust.jcommander.Parameters;

import java.util.ArrayList;
import java.util.List;

/**
 * Arguments or randomizer2.
 */
@Parameters(commandDescription = "Randomize one or several sbi/sbip files.")

public class Randomizer2Arguments {
    @Parameter(required = true, variableArity = true, names = {"-i", "--input-files"}, description = "Input files in .bsi/.bsip format.")
    public List<String> inputFiles = new ArrayList<>();

    @Parameter(required = false, names = {"-b", "--records-per-bucket"}, description = "Number of records to store in each bucket.")
    int recordsPerBucket = 20000;

    @Parameter(required = true, names = {"-o", "--output-prefix"}, description = "Prefix for the output filenames.")
    public String outputFile;

    @Parameter(required = false, names = {"-c", "--chunk-size"}, description = "Size of chunks for each bucket writer.")

    public int chunkSizePerWriter = 1000;
}

