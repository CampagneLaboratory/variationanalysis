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
    @Parameter( names = { "--score-indels"}, description = "When this flag is provided, we score indel alleles as" +
            " we would SNPs. When this flag is omitted, we trim alleles longer than 1 base to exactly one base before" +
            "comparing for allele match. For instance, with --score-indel, we consider A--------/ATTATTTAT vs  A/A to "+
            " be a mismatch. Without it, we trim the indel to A/A and find a match.")
    boolean scoreIndels;

}