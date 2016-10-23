package org.campagnelab.dl.varanalysis.tools;

import com.beust.jcommander.JCommander;
import com.beust.jcommander.ParameterException;

/**
 * Created by fac2003 on 10/23/16.
 */
public abstract class AbstractTool {
    ToolArguments arguments;

    protected void parseArguments(String[] args, String commandName, ToolArguments arguments) {
        this.arguments = arguments;
        JCommander commander = new JCommander(arguments);
        commander.setProgramName(commandName);
        try {
            commander.parse(args);
        } catch (ParameterException e) {

            commander.usage();
            System.out.flush();
            throw e;
        }
    }

    public abstract void execute();
}
