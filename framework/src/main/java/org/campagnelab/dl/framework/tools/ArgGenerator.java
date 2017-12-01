package org.campagnelab.dl.framework.tools;

import com.beust.jcommander.JCommander;
import org.apache.commons.io.FileUtils;
import org.campagnelab.dl.framework.tools.arguments.ArgGeneratorArguments;

import java.io.*;
import java.util.*;

/**
 * Created by rct66 on 1/5/17.
 * Tool to generate random argument combinations.
 */
public class ArgGenerator {

    //this maps argument names to a list of options or (or range boundaries) for the argument.
    LinkedHashMap<String, List<String>> options = new LinkedHashMap();
    //this maps argument names to a type: either categorical, int, or float.
    Map<String, String> types = new HashMap();


    public static void main(String[] args) {


        ArgGeneratorArguments arguments = new ArgGeneratorArguments();
        JCommander commander = new JCommander(arguments);
        commander.setProgramName("ArgGenerator");
        commander.parse(args);


        ArgGenerator driver = new ArgGenerator(new Date().getTime());
        try {
            driver.configure(arguments.argConfig);
            FileUtils.writeStringToFile(new File(arguments.outputFilename), driver.generateCommands(arguments.numCommands));

        } catch (IOException e) {
            e.printStackTrace();
            throw new RuntimeException("There was a problem parsing the configuration file.");
        }

    }

    private Random rand;

    ArgGenerator(long seed) {
        rand = new Random(seed);
    }


    protected void configure(String configPath) throws IOException {
        String configString = FileUtils.readFileToString(new File(configPath));
        configureWithString(configString);
    }

    protected void configureWithString(String configString) throws IOException {

        String argName;
        BufferedReader br = new BufferedReader(new StringReader(configString));
        while ((argName = br.readLine()) != null) {
            List<String> options = new ArrayList<String>();
            //store type in a map
            types.put(argName, br.readLine());

            //now iterate over arg options to make a list and put it in the option map
            String argOption;
            while ((argOption = br.readLine()) != null) {
                if (argOption.equals("")) {
                    break;
                }
                options.add(argOption);
            }
            this.options.put(argName, options);
        }
        br.close();
    }

    protected String generateCommands(int numCommands) throws IOException {

        StringBuffer result = new StringBuffer();

        for (int i = 0; i < numCommands; i++) {
            StringBuffer command = new StringBuffer();
            for (Map.Entry<String, List<String>> entry : options.entrySet()) {
                String argName = entry.getKey();
                String argType = types.get(argName);
                List<String> argOptions = entry.getValue();
                String format = "%s %s";
                String value;
                switch (argType) {
                    case "categorical":
                        value = argOptions.get(rand.nextInt(argOptions.size()));
                        break;
                    case "uniform":
                        float minF = Float.parseFloat(argOptions.get(0));
                        float maxF = Float.parseFloat(argOptions.get(1));
                        value = Float.toString((rand.nextFloat() * (maxF - minF) + minF));
                        break;
                    case "log-uniform":
                        double minD = Double.parseDouble(argOptions.get(0));
                        double maxD = Double.parseDouble(argOptions.get(1));
                        if (minD == 0 || maxD == 0) {
                            System.err.printf("%s: you cannot use a zero bound with log-uniform.", argName);
                            System.exit(1);
                        }
                        double minLog = Math.log(minD);
                        double maxLog = Math.log(maxD);
                        double valueLog = (rand.nextDouble() * (maxLog - minLog) + minLog);
                        double valueD = Math.exp(valueLog);
                        value = Double.toString(valueD);
                        break;
                    case "flag":
                        // decide if we should include this flag argument:
                        if (rand.nextBoolean()) {
                            format = "%s";
                            value = "";
                        } else {
                            argName = "";
                            format = "";
                            value = "";
                        }
                        break;
                    case "int":
                        int minI = Integer.parseInt(argOptions.get(0));
                        int maxI = Integer.parseInt(argOptions.get(1));
                        value = Integer.toString((rand.nextInt(maxI - minI + 1) + minI));
                        break;
                    default:
                        throw new RuntimeException("There was a problem parsing the config. A non-existent argType ([categorical|uniform|log-uniform|int]) may have been used: " + argType);
                }
                command.append(String.format(format, argName, value));
                command.append(" ");
            }
            //remove trailing space
            result.append(command.toString().trim());


            if (i!=numCommands-1){
                result.append("\n");
            }
        }
        return result.toString();
    }


}
