package org.campagnelab.dl.genotype.tools;

import com.beust.jcommander.Parameter;
import org.campagnelab.dl.framework.tools.arguments.ToolArguments;
import org.campagnelab.dl.genotype.mappers.GenotypeMapperV37;

import java.util.ArrayList;
import java.util.List;

public class WriteVectorsArguments implements ToolArguments {
    @Parameter(required = true, names = {"-t", "--training-set"},
            description = "Training set, must be provided in .sbi/.sbip format (produced with Goby3).")
    public String trainingSet = null;

    @Parameter(names = "--feature-mapper",
            description = "Fully qualified name of the feature mapper class.")
    public String featureMapperClassname = GenotypeMapperV37.class.getName();
//
//    @Parameter(names = "--label-mapper", required = true,
//            description = "Fully qualified name of the label mapper class.")
//    public String labelMapperClassname = null;

    @Parameter(names = {"-o","--output"}, required = true,
            description = "Output basename (possibly including path where output files will be created (i.e., /a/b/out).")
    public String outputBasename;
}
