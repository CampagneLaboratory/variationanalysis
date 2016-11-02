package org.campagnelab.dl.varanalysis.tools;

import com.beust.jcommander.Parameter;
import org.campagnelab.dl.model.utils.mappers.FeatureMapperV18;

import java.util.ArrayList;
import java.util.List;

/**
 * Arguments for MapFeatures.
 * Created by fac2003 on 11/2/16.
 */
public class MapFeaturesArguments implements ToolArguments {
    @Parameter(required = true, names = {"-t", "--training-sets"}, variableArity = true, description = "Training sets, must be provided in .sbi/.sbip format (produced with Goby3). When more than one dataset is provided (multiple -t options), the " +
            "datasets are concatenated.")
    public List<String> trainingSets = new ArrayList<>();

    @Parameter(names = "--feature-mapper", description = "Fully qualified name of the feature mapper class.")
    public String featureMapperClassname = FeatureMapperV18.class.getCanonicalName();

    @Parameter(names = "--mini-batch-size", description = "The size of the training minibatch")
    public int miniBatchSize = 32;

    @Parameter(names = "--trio", description = "Use to train trio models. The training and validation datasets must have three samples, parents first, patient last.")
    public boolean isTrio = false;

    public String[] getTrainingSets() {
        return this.trainingSets.toArray(new String[this.trainingSets.size()]);
    }

    @Parameter(names = {"-o","--output"},required = true,description = "Output basename (possibly including path where output files will be created (i.e., /a/b/out).")
    public
    String outputBasename;

    @Parameter(names = {"-n","--write-n"}, description = "Write at most n records, then stop.")
    public long writeAtMostN = Long.MAX_VALUE;
}
