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
    @Parameter(names = {"--het-filter"}, description = "Keep sites that are either HET, HOM or both (ALL).")
    Show showFilter = Show.ALL;

    enum Show {
        HET, HOM, ALL
    }

    @Parameter(names = {"--format"}, description = "Report format (TABLE or VCF).")
    OutputFormat outputFormat = OutputFormat.TABLE;

    enum OutputFormat {
        TABLE, VCF
    }

    @Parameter(names = {"--only-variants"}, description = "Only show sites where one allele is not reference.")
    boolean onlyVariants;

    @Parameter(names = {"--num-variants-expected"}, description = "The number of variants expected in the test set, " +
            "used to estimate the false negative rate. Can be an estimate since we use max(num-variants-expected,TP) to " +
            "calculate TP+FN.")
    int numVariantsExpected;

    @Parameter(names = {"--minimum-coverage"}, description = "The minimum coverage needed to report a site, " +
            "used to filter exome results where some off-target hits are expected with very low coverage (e.g., 10). Sites with" +
            "at least the number of reads mapping are reported. Default 0 (no coverage filter)" )
    int minimumCoverage=0;
}