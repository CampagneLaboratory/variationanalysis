package org.campagnelab.dl.genotype.tools;

import com.beust.jcommander.Parameter;
import org.campagnelab.dl.framework.tools.arguments.ToolArguments;

import java.util.ArrayList;
import java.util.List;

/**
 * Arguments for the {@link SBISimulator} tool.
 */
public class SBISimulatorArguments implements ToolArguments {

    @Parameter(required = true, names = {"-o", "--output-filename"}, description = "Output filename for the sbi.")
    public String outputFilename;

    @Parameter(required = true, names = {"-i", "--input-file"}, variableArity = true,
            description = "Input file with the variant map, must be provided in .varmap format (produced with Goby3).")
    public String inputFile = null;

    @Parameter(required = false, names = {"--include"}, variableArity = true,
            description = "Chromosome(s) to include.")
    public List<String> includeChromosomes = new ArrayList<>();

    @Parameter(required = false, names = {"--exclude"}, variableArity = true,
            description = "Chromosome(s) to exclude.")
    public List<String> excludeChromosomes = new ArrayList<>();


    @Parameter(names = "--read-N", description = "Read at most N chromosomes from the varmap, then stop.")
    public int readN = Integer.MAX_VALUE;

    @Parameter(names = "--genomic-context-length", description = "Length of the genomic context (odd number, 1,3+).")
    public int genomicContextLength=1;

    @Parameter(required = true, names = {"--genome"}, description = "Basename of a goby indexed genome.")
    public String genome;

    @Parameter(names ={"-v", "--verbose"}, description = "Be more verbose.")
    boolean verbose;
}

