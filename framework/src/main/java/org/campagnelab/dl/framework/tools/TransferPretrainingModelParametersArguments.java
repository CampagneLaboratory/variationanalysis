package org.campagnelab.dl.framework.tools;

import com.beust.jcommander.Parameter;
import org.campagnelab.dl.framework.tools.arguments.RecordingToolArguments;
import org.campagnelab.dl.framework.tools.arguments.ToolArguments;

/**
 * Created by joshuacohen on 12/19/16.
 */
public abstract class TransferPretrainingModelParametersArguments extends TrainingArguments {
    @Parameter(names = "--pretraining-model-path", required = true, description = "Path to pretrained model to load to initialize parameters.  ")
    public String pretrainingModelPath;
    @Parameter(names = "--pretraining-model-name", required = true, description = "Name of pretrained model to load to initialize parameters. ")
    public String pretrainingModelName;
    @Parameter(names = "--model-path", description = "Directory to save created model in. ")
    public String modelPath = null;
    @Parameter(names = "--model-prefix", description = "Prefix for saved model. ")
    public String modelPrefix = null;

    public abstract String domainDescriptorName();
}
