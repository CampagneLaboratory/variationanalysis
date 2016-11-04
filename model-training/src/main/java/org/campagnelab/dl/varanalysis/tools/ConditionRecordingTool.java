package org.campagnelab.dl.varanalysis.tools;

import com.beust.jcommander.JCommander;
import com.beust.jcommander.ParameterDescription;
import com.beust.jcommander.ParameterException;

import java.io.IOException;
import java.util.HashMap;
import java.util.Map;

/**
 * A tool that records model conditions and results to a file.
 */
public abstract class ConditionRecordingTool<T extends ToolArguments> extends AbstractTool<T> {
    /**
     * Stores statistics produced by this tool.
     */
    private Map<String, Object> resultValues;

    /**
     * Get the map of result names to values
     * @return map of results to values
     */
    public Map<String, Object> resultValues() {
        return resultValues;
    }
    // Map of command-line strings for fields to their set values if they were specified on the command-line
    private Map<String, Object> setFieldValues;
    // Map of command-line strings for fields to their default values if they weren't specified
    private Map<String, Object> defaultFieldValues;

    /**
     * Store the command-line strings for fields with their values after argument parsing
     *
     * @param commander JCommander instance with parsed arguments
     */
    public void storeFieldValues(JCommander commander) {
        setFieldValues = new HashMap<>();
        defaultFieldValues = new HashMap<>();
        for (ParameterDescription p : commander.getParameters()) {
            if (p.isAssigned()) {
                setFieldValues.put(p.getLongestName(), p.getParameterized().get(this));
            } else {
                defaultFieldValues.put(p.getLongestName(), p.getParameterized().get(this));
            }
        }
    }

    /**
     * Get the command-line strings for fields with their set values if they were specified on the command line
     *
     * @return Map from command-line strings for fields to values through command-line arguments
     */
    public Map<String, Object> setFieldValues() {
        return setFieldValues;
    }

    /**
     * Get the command-line strings for fields with their set values if they were set by default values
     *
     * @return Map from command-line strings for fields to values by default
     */
    public Map<String, Object> defaultFieldValues() {
        return defaultFieldValues;
    }

    public void writeModelingConditions(RecordingToolArguments arguments) throws IOException {
        String header = "Results|Specified_Arguments|Default_Arguments\n";
        ModelConditionHelper.createLogFile(arguments.modelConditionFilename, header);
        ModelConditionHelper.appendToLogFile(arguments.modelConditionFilename,
                ModelConditionHelper.fieldMapToString(resultValues()),
                ModelConditionHelper.fieldMapToString(setFieldValues()),
                ModelConditionHelper.fieldMapToString(defaultFieldValues()));

    }
    protected void parseArguments(String[] args, String commandName, T arguments) {
        this.arguments = arguments;
        JCommander commander = new JCommander(arguments);
        commander.setProgramName(commandName);
        try {
            commander.parse(args);
            storeFieldValues(commander);
        } catch (ParameterException e) {

            commander.usage();
            System.out.flush();
            throw e;
        }
    }

}
