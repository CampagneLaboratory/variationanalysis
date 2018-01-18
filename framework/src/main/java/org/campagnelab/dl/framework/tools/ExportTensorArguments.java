package org.campagnelab.dl.framework.tools;

import com.beust.jcommander.Parameter;
import org.campagnelab.dl.framework.tools.arguments.ToolArguments;

import java.util.ArrayList;
import java.util.List;
import java.util.Set;

/**
 * Arguments for ExportTensor.
 * Created by fac2003 on 1/16/2018.
 */
public class ExportTensorArguments implements ToolArguments {
    @Parameter(required = true, names = {"-t", "--training-sets"}, variableArity = true, description = "Training sets, must be provided in .sbi/.sbip format (produced with Goby3). When more than one dataset is provided (multiple -t options), the " +
            "datasets are concatenated.")
    public List<String> trainingSets = new ArrayList<>();

    @Parameter(names = "--feature-mapper", required = true, description = "Fully qualified name of the feature mapper class.")
    public String featureMapperClassname = null;

    @Parameter(names = "--label-mapper", required = false, description = "Fully qualified name of the label mapper class.")
    public String labelMapperClassname = null;

    @Parameter(names = "--mini-batch-size", description = "The size of the training minibatch")
    public int miniBatchSize = 32;

    @Parameter(names = {"n","--export-n"}, description = "Export at most n records.")
    public int exportN=Integer.MAX_VALUE;

    public String[] getTrainingSets() {
        return this.trainingSets.toArray(new String[this.trainingSets.size()]);
    }

    @Parameter(names = {"-o","--output"},required = false, description = "Output basename (possibly including path where output files will be created (i.e., /a/b/out).")
    public
    String outputBasename="exported-tensors";

    @Parameter(names = {"--export-inputs"}, description = "Restrict the exports to the inputs with the specified names. " +
            "If you don't know the names of the inputs, run the tool without this argument, and all input names will be shown.")
    public Set<String> inputNamesToExport;

    @Parameter(names = {"--export-outputs"}, description = "Restrict the exports to the inputs with the specified names. " +
            "If you don't know the names of the inputs, run the tool without this argument, and all input names will be shown.")
    public Set<String> outputNamesToExport;

    @Parameter(required = true, names = "--sample-names",
            description = "Names for each of the samples in the reads", variableArity = true)
    public List<String> sampleNames = new ArrayList<>();

    @Parameter(required = true, names = "--sample-types",
            description = "Types for each of the samples in the reads", variableArity = true)
    public List<String> sampleTypes = new ArrayList<>();

    @Parameter(required = true, names = "--sample-ids",
            description = "IDs for each of the samples in the reads", variableArity = true)
    public List<Integer> sampleIds = new ArrayList<>();

}
