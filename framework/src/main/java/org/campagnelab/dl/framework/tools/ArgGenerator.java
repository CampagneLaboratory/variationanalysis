package org.campagnelab.dl.framework.tools;

import com.beust.jcommander.JCommander;
import org.campagnelab.dl.framework.tools.arguments.ArgGeneratorArguments;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintWriter;
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

    Random rand = new Random();


    public static void main(String[] args) {


        ArgGeneratorArguments arguments = new ArgGeneratorArguments();
        JCommander commander = new JCommander(arguments);
        commander.setProgramName("ArgGenerator");
        commander.parse(args);


        ArgGenerator driver = new ArgGenerator();
        try {
            driver.configure(arguments.argConfig);
            driver.generateCommands(arguments.outputFilename, arguments.numCommands);
        } catch (IOException e) {
            e.printStackTrace();
            throw new RuntimeException("There was a problem parsing the configuration file.");
        }

    }

    ArgGenerator() {
    }

    ;


    void configure(String configPath) throws IOException {

        BufferedReader br = new BufferedReader(new FileReader(configPath));
        String argName;
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

    void generateCommands(String outputPath, int numCommands) throws IOException {

        PrintWriter writer = new PrintWriter(outputPath, "UTF-8");

        for (int i = 0; i < numCommands; i++) {
            StringBuffer command = new StringBuffer();
            for (Map.Entry<String, List<String>> entry : options.entrySet()) {
                String argName = entry.getKey();
                String argType = types.get(argName);
                List<String> argOptions = entry.getValue();

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
                    case "int":
                        int minI = Integer.parseInt(argOptions.get(0));
                        int maxI = Integer.parseInt(argOptions.get(1));
                        value = Integer.toString((rand.nextInt(maxI - minI + 1) + minI));
                        break;
                    default:
                        throw new RuntimeException("There was a problem parsing the config. A non-existent argType ([categorical|uniform|log-uniform|int]) may have been used.");
                }
                command.append(argName + " " + value + " ");
            }
            //remove trailing space
            writer.println(command.toString().trim());
        }
        writer.close();
    }


}
