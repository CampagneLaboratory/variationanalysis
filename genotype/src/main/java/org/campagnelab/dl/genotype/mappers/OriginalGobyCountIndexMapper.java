package org.campagnelab.dl.genotype.mappers;

import org.campagnelab.dl.framework.mappers.FeatureNameMapper;
import org.campagnelab.dl.framework.mappers.OneHotHashModuloMapper;
import org.campagnelab.dl.varanalysis.protobuf.BaseInformationRecords;

/**
 * Map the original index of the goby count in the sbi.
 * Created by fac2003 on 2/19/17.
 */
public class OriginalGobyCountIndexMapper extends OneHotHashModuloMapper<BaseInformationRecords.BaseInformationOrBuilder> implements FeatureNameMapper<BaseInformationRecords.BaseInformationOrBuilder> {
    int genotypeIndex;
    int sampleIndex;

    public OriginalGobyCountIndexMapper(int sampleIndex, int genotypeIndex) {
        super(20, record -> record.getSamples(sampleIndex).getCounts(genotypeIndex).getGobyGenotypeIndex());
    }

    @Override
    public String getFeatureName(int featureIndex) {
        return "originalGobyCountInex"+featureIndex;
    }
}
