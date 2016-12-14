package org.campagnelab.dl.genotype.tools;

import org.campagnelab.dl.framework.domains.DomainDescriptor;
import org.campagnelab.dl.framework.tools.Show;
import org.campagnelab.dl.genotype.learning.domains.GenotypeDomainDescriptor;
import org.campagnelab.dl.somatic.tools.SomaticShowArguments;
import org.campagnelab.dl.varanalysis.protobuf.BaseInformationRecords;

import java.util.function.Function;

/**
 * Show somatic/sbi records.
 */
public class ShowG extends Show<BaseInformationRecords.BaseInformation, GenotypeShowArguments> {
    public static void main(String[] args) {

        Show show = new ShowG();
        show.parseArguments(args, "ShowG", show.createArguments());
        show.execute();
    }

    @Override
    protected DomainDescriptor<BaseInformationRecords.BaseInformation> domainDescriptor() {
        return new GenotypeDomainDescriptor(args().modelPath);
    }

    @Override
    public GenotypeShowArguments createArguments() {
        return new GenotypeShowArguments();
    }

    @Override
    protected Function<BaseInformationRecords.BaseInformation, String> getConverter(String reportType) {
        return args().getConverter();
    }
}
