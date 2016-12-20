package org.campagnelab.dl.genotype.tools;

import com.beust.jcommander.Parameter;
import com.beust.jcommander.Parameters;
import org.campagnelab.dl.framework.tools.arguments.RecordingToolArguments;

import java.util.ArrayList;
import java.util.Date;
import java.util.List;

/**
 * Arguments for tools that train models.
 * Created by fac2003 on 8/20/16.
 */
@Parameters(commandDescription = "Print the feature names of a specific feature mapper")

public class FeatureNameArguments extends RecordingToolArguments {


    @Parameter(names = "--feature-mapper", description = "Fully qualified name of the feature mapper class.")
    public String featureMapperClassname = "org.campagnelab.dl.genotype.mappers.GenotypeMapperV4";


}
