package org.campagnelab.dl.somatic.tools;

import com.beust.jcommander.Parameter;
import com.beust.jcommander.Parameters;
import org.campagnelab.dl.framework.tools.arguments.ToolArguments;

import java.util.List;

/**
 * Created by fac2003 on 9/2/16.
 */
@Parameters(commandDescription = "Split a file into several components.")

public class SplitArguments implements ToolArguments {
    @Parameter(required = true, names = {"-i", "--input-file"}, description = "Input file in .bsi/.bsip format.")
    public String inputFile;

    @Parameter(required = true, names = {"-f", "--fraction"}, description = "Fraction of the input file to put in a destination Fractions will be normalized before use so you" +
            "may use -f 0.3 and -f 0.7 or -f 70 and -f 30.")
    List<Float> fractions;

    @Parameter(required=true, names = {"-s", "--suffix"}, description = "Suffix that will be added to filename for this fraction.")
    List<String> suffixes;

    @Parameter(required=true, names = {"-o", "--output-prefix"}, description = "Prefix for the output filenames.")
    public String outputFile;

    @Parameter(required=false, names = {"-n", "--write-n"}, description = "Write at most n records.")
    public long writeN=Long.MAX_VALUE;

    @Parameter(required=false, names = { "--destination-override"}, description = "When provided, a syntax such as test:chr1,chr2 will force all records with chromosome chr1 or chr2 to be written to the test suffix output.")
    public String destinationOverride=null;

    @Parameter(required=false, names = { "--random-seed"}, description = "Seed for random generator used to determine in which file to output a site.")
    long randomSeed=232323;

}

