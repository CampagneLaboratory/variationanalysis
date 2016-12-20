package org.campagnelab.dl.framework.tools;

import com.beust.jcommander.Parameter;
import org.campagnelab.dl.framework.tools.arguments.ToolArguments;

/**
 * Created by joshuacohen on 12/19/16.
 */
public class TransferPretrainingModelParametersArguments implements ToolArguments {
    @Parameter(names = "--pretraining-model-path", description = "If provided, use the pretraining model at pretrainedModelPath to initialize parameters. ")
    public String pretrainingModelPath = null;
    @Parameter(names = "--pretraining-model-name", description = "The name of the pretrained model to load to initialize parameters. ")
    public String pretrainingModelName = "best";
    @Parameter(names = "--eos-character", description = "If provided, use as EOS character index for alignment. If not, adds EOS to the input vocab. Must be specified if specified during pretraining, and likewise if not. " )
    public Integer eosIndex = null;
    @Parameter(names = "--model-path", required = true, description = "Directory to save created model in. ")
    public String modelPath;
    @Parameter(names = "--model-prefix", required = true, description = "Prefix for saved model. ")
    public String modelPrefix;
}
