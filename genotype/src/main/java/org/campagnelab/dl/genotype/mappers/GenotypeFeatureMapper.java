package org.campagnelab.dl.genotype.mappers;

import org.campagnelab.dl.framework.mappers.ConfigurableFeatureMapper;
import org.campagnelab.dl.somatic.mappers.NamingConcatFeatureMapper;
import org.campagnelab.dl.varanalysis.protobuf.BaseInformationRecords;

/**
 * Created by fac2003 on 12/16/16.
 */
public abstract class GenotypeFeatureMapper extends NamingConcatFeatureMapper<BaseInformationRecords.BaseInformationOrBuilder> implements ConfigurableFeatureMapper {
    public boolean sortCounts;
    public boolean withDistinctAlleleCounts;
    public boolean withCombinedLayer;
    public boolean hasIsVariantLabelMapper;
    protected int genomicContextLength;

    public void setGenomicContextLength(int genomicContextLength) {
        if ((genomicContextLength % 2)==0) {
            throw new RuntimeException("The genomic context length must be an odd number, usually in the range 21-61");
        }
        this.genomicContextLength = genomicContextLength;
    }
}
