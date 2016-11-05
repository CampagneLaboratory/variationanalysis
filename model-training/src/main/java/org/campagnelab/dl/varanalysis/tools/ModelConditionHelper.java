package org.campagnelab.dl.varanalysis.tools;

import org.apache.commons.io.FileUtils;
import org.apache.commons.lang3.StringUtils;

import java.io.*;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;

/**
 * Methods to help write model condition files.
 * Created by joshuacohen on 11/2/16.
 */
public class ModelConditionHelper {
    /**
     * Output string for field maps in format that can be accepted as command line arguments
     * @param fieldMap map containing fields as command-line strings and values
     * @return String representation of this map
     */
    public static String fieldMapToString(Map<String, Object> fieldMap) {
       assert fieldMap!=null: "fieldMap must not be null";
        List<String> fieldMapStrings = new LinkedList<>();
        for (Map.Entry<String, Object> fieldValue: fieldMap.entrySet()) {
            fieldMapStrings.add(fieldValue.getKey() + " " + fieldValue.getValue());
        }
        return StringUtils.join(fieldMapStrings, ' ');
    }

    /**
     * Set up log file for use
     * @param logPath full path to log file
     * @param header header for log file
     * @throws IOException
     */
    public static void createLogFile(String logPath, String header) throws IOException {
        File logFile = new File(logPath);
        if (!logFile.exists()) {
            String logPathParent = logPath.substring(0, logPath.lastIndexOf("/"));
            FileUtils.forceMkdir(new File(logPathParent));
            Writer logFileWriter = new BufferedWriter(new FileWriter(logFile));
            logFileWriter.append(header);
            logFileWriter.close();
        }

    }

    /**
     * Append strings to the given log file
     * @param logPath full path to log file
     * @param lineArgs strings to add to the log file
     * @throws IOException
     */
    public static void appendToLogFile(String logPath, String tag, String... lineArgs) throws IOException {
        Writer logFileWriter = new BufferedWriter(new FileWriter(logPath, true));
        logFileWriter.append(tag);
        logFileWriter.append('|');
        logFileWriter.append(StringUtils.join(lineArgs, '|') + "\n");
        logFileWriter.close();
    }
}
