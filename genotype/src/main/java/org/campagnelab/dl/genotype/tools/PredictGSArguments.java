package org.campagnelab.dl.genotype.tools;

import com.beust.jcommander.Parameter;

/**
 * Arguments for the {@link PredictGS} tool.
 *
 * @author manuele
 */
public class PredictGSArguments extends PredictGArguments {

    @Parameter(names = {"--split-indels"}, description = "Write indel sites in a separate VCF.")
    boolean splitIndels = false;


}
