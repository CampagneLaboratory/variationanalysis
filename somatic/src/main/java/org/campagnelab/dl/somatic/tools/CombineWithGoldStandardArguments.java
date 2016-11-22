package org.campagnelab.dl.somatic.tools;

import com.beust.jcommander.Parameter;
import org.campagnelab.dl.framework.tools.arguments.AbstractTool;
import org.campagnelab.dl.framework.tools.arguments.ToolArguments;

import java.util.List;

/**
 * Created by fac2003 on 11/22/16.
 */
public class CombineWithGoldStandardArguments implements ToolArguments {
    @Parameter(required = true, names = {"-i", "--input-file"}, description = "Input file in .bsi/.bsip format.")
    String sbiFilename;

    @Parameter(required = true, names = {"-a", "--annotations"}, description = "Annotations in TSV format: chromosome\tposition.")
    String annotationFilename;
    @Parameter(required = true, names = {"-o", "--output-file"}, description = "Output file in .bsi/.bsip format.")
    public String outputFilename;

    @Parameter(required = true, names = {"-f", "--sampling-fraction"}, description = "Fraction of the input file to write to the output. Annotated site are written irrespective of fraction.")
    float samplingFraction=1f;
}
