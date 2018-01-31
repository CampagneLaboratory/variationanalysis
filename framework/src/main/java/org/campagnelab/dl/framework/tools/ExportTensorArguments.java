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
    @Parameter(required = true, names = {"-i",}, variableArity = true, description = "Input files sets, must be provided " +
            "in .sbi/.sbip format (produced with Goby3). " +
            "When more than one dataset is provided (multiple -i options), the " +
            "datasets are concatenated.")
    public List<String> trainingSets = new ArrayList<>();

    @Parameter(names = "--feature-mapper", required = true, description = "Fully qualified name of the feature mapper class.")
    public String featureMapperClassname = null;

    @Parameter(names = "--mini-batch-size", description = "The size of the minibatch used to map features. Larger sizes can be more efficient but require more memory.")
    public int miniBatchSize = 64;

    @Parameter(names = {"-n", "--export-n"}, description = "Export at most n records.")
    public int exportN = Integer.MAX_VALUE;

    public String[] getTrainingSets() {
        return this.trainingSets.toArray(new String[this.trainingSets.size()]);
    }

    @Parameter(names = {"-o", "--output"}, required = false, description = "Output basename (possibly including path where output files will be created (i.e., /a/b/out).")
    public
    String outputBasename = "./exported-tensors";

    @Parameter(names = {"--export-input"}, description = "Restrict the exports to the inputs with the specified names. " +
            "If you don't know the names of the inputs, run the tool without this argument, and all input names will be shown. Repeat the option " +
            "as necessary with different input names.")
    public Set<String> inputNamesToExport;

    @Parameter(names = {"--export-output"}, description = "Restrict the exports to the inputs with the specified names. " +
            "If you don't know the names of the inputs, run the tool without this argument, and all input names will be shown. Repeat the option " +
            "as necessary with different output names.")
    public Set<String> outputNamesToExport;

    @Parameter(required = true, names = "--sample-name", description = "Name of the samples in the input file(s). " +
            "One name per sample, repeat the option as necessary.",
            variableArity = true)
    public List<String> sampleNames = new ArrayList<>();

    @Parameter(required = true, names = "--sample-type",
            description = "Type of the samples in the input file(s). " +
                    "The type is used to represent a property about the sample, such as germline or somatic. " +
                    "One type per sample, repeat the option as necessary.", variableArity = true)
    public List<String> sampleTypes = new ArrayList<>();

    @Parameter(required = false, names = "--sample-index",
            description = "Index of a sample to export. If an .sbi file is created with several samples, you can choose which sample to export with this option. Repeat the option as necessary with different indices." +
                    "Defaults to 0 to export the first sample only.", variableArity = true)
    public List<Integer> sampleIds = defaultSampleIndices();

    @Parameter(names = "--ploidy", description = "The organism ploidy (2 for humans, more for some plants). This parameter controls some of the mapped features and dimension of" +
            "outputs. ")
    public int ploidy = 2;

    @Parameter(names = { "--genomic-context-length"}, description = "Length of genomic context to use around site, in mapped features. The larger the number," +
            "the more bases are presented to the neural net around the genomic site of interest, the more features are used to represent the genomic context." +
            "This parameter controls some of the mapped features and the dimension of the input feature vector(s).")
    public int genomicContextLength=29;

    @Parameter(names = "--label-smoothing-epsilon", description = "Value of epsilon for label smoothing. Zero (default) is no smoothing. Try small values (<0.1)." +
            "This parameter controls some of the mapped outputs. ")
    public float labelSmoothingEpsilon = 0;

    private List<Integer> defaultSampleIndices() {
        List<Integer> result = new ArrayList<>();
        result.add(0);
        return result;
    }

    @Parameter(names = "--vector-file-type", description = "Type of .vec file to write out. Can be 'text' or 'binary' vec format. " +
            "Text format is compressed with gzip and smaller for storage and transport. Binary format is useful for random access and is used for " +
            "caching data for faster training with shuffled datasets.")
    public String vecFileType = "text";

}
