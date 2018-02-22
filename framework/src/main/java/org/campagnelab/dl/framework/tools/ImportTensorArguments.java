package org.campagnelab.dl.framework.tools;

import com.beust.jcommander.Parameter;
import org.campagnelab.dl.framework.tools.arguments.ToolArguments;

import java.util.List;
import java.util.function.Consumer;

public class ImportTensorArguments implements ToolArguments{
    @Parameter(required = true, names = {"-i", "--input"}, description = "Input .vec file")
    public String inputPath;

    @Parameter(names = {"-n", "--import-n"}, description = "Import at most n records.")
    public long importN = Long.MAX_VALUE;

    @Parameter(names = {"-s", "--sample-id"}, description = "Sample ID to read in records for")
    public int sampleId = 0;

    @Parameter(required = true, names = {"-v", "--vector-names"}, variableArity = true,
            description = "Vector names to read in vectors for")
    public List<String> vectorNames;

    @Parameter(names = "--mini-batch-size", description = "The size of the minibatch")
    public int miniBatchSize = 1;

    public Consumer<VectorReader.RecordVectors> processVectors() {
        return System.out::println;
    }
}
