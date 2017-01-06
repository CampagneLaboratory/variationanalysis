package org.campagnelab.dl.framework.tools.arguments;
/**
 * A version of ToolArguments that records command line parameters and results, then writes both to a log file.
 */

import com.beust.jcommander.Parameter;

public class ArgGeneratorArguments {
    @Parameter( names = {"--config"}, description = "Path to the argument configuration file.", required = true)
    public String argConfig;

    @Parameter( names = {"--output"}, description = "output file which will contain generated commands, each in new line.", required = true)
    public String outputFilename;

    @Parameter( names = {"--num-commands"}, description = "number of commands to generate and write to output file.")
    public int numCommands = 10;

}