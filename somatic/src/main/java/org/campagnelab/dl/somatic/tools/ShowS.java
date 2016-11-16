package org.campagnelab.dl.somatic.tools;

import org.campagnelab.dl.framework.domains.DomainDescriptor;
import org.campagnelab.dl.framework.tools.Show;
import org.campagnelab.dl.somatic.learning.domains.SomaticMutationDomainDescriptor;
import org.campagnelab.dl.varanalysis.protobuf.BaseInformationRecords;

import java.util.function.Function;

/**
 * Show somatic/sbi records.
 */
public class ShowS extends Show<BaseInformationRecords.BaseInformation, SomaticShowArguments> {
    public static void main(String[] args) {

        Show show = new ShowS();
        show.parseArguments(args, "ShowS", show.createArguments());
        show.execute();
    }

    @Override
    protected DomainDescriptor<BaseInformationRecords.BaseInformation> domainDescriptor() {
        return new SomaticMutationDomainDescriptor(args().modelPath);
    }

    @Override
    public SomaticShowArguments createArguments() {
        return new SomaticShowArguments();
    }

    @Override
    protected Function<BaseInformationRecords.BaseInformation, String> getConverter(String reportType) {
        return args().getConverter();
    }
}
