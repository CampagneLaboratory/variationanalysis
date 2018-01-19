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
    public boolean withCombinedLayerRef;
    public boolean hasIsVariantLabelMapper;
    public boolean withSoftmaxGenotype;
    public static int MAX_GENOTYPES = 3;

    protected int sampleIndex;
    /**
     * Set the sampleIndex before calling configure. Sample index cannot be changed after configure has been called.
     * @param sampleIndex Index of the sample in the sbi record to use for mapping features.
     */
    public void setSampleIndex(int sampleIndex) {
        this.sampleIndex = sampleIndex;
    }

}
