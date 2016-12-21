package org.campagnelab.dl.genotype.tools;


import com.beust.jcommander.Parameter;
import org.campagnelab.dl.framework.tools.PredictArguments;

/**
 * Example of Predict implementation. This class performs predictions with a model trained by TrainModelS.
 *
 * @author Remi Torracinta
 *         Created by rct66 on 12/7/16.
 */
public class PredictGArguments extends PredictArguments {
    @Parameter( names = { "--het-filter"}, description = "Keep sites that are either HET, HOM or both (ALL).")
    Show showFilter=Show.ALL;
    enum Show{
        HET, HOM, ALL
    }

}