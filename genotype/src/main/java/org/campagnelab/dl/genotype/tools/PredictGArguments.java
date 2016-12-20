package org.campagnelab.dl.genotype.tools;


import com.beust.jcommander.Parameter;
import org.campagnelab.dl.framework.domains.prediction.Prediction;
import org.campagnelab.dl.framework.tools.Predict;
import org.campagnelab.dl.framework.tools.PredictArguments;
import org.campagnelab.dl.genotype.learning.domains.predictions.HomozygousPrediction;
import org.campagnelab.dl.genotype.learning.domains.predictions.SingleGenotypePrediction;
import org.campagnelab.dl.genotype.predictions.GenotypePrediction;
import org.campagnelab.dl.varanalysis.protobuf.BaseInformationRecords;

import java.io.PrintWriter;
import java.util.List;

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