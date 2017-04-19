package org.campagnelab.dl.genotype.mappers;

import org.campagnelab.dl.framework.mappers.FeatureMapper;
import org.campagnelab.dl.framework.mappers.NoMaskFeatureMapper;
import org.campagnelab.dl.somatic.mappers.NoMasksLabelMapper;
import org.campagnelab.dl.varanalysis.protobuf.BaseInformationRecords;
import org.nd4j.linalg.api.ndarray.INDArray;

/**
 * A label mapper that can optionally pre-sort genotypes by decreasing count.
 * Created by fac2003 on 12/20/16.
 *
 */
public abstract class CountSortingLabelMapper extends NoMasksLabelMapper<BaseInformationRecords.BaseInformation> {
    int[] indices = new int[]{0, 0};
    private final boolean sortCounts;

    public CountSortingLabelMapper(boolean sortCounts) {
        this.sortCounts = sortCounts;
    }

    @Override
    public void mapLabels(BaseInformationRecords.BaseInformation record, INDArray labels, int indexOfRecord) {
        indices[0] = indexOfRecord;

        for (int labelIndex = 0; labelIndex < numberOfLabels(); labelIndex++) {
            indices[1] = labelIndex;
            labels.putScalar(indices, produceLabel(sortedCountRecord, labelIndex));
        }
    }

    @Override
    public void prepareToNormalize(BaseInformationRecords.BaseInformation record, int indexOfRecord) {
        if (sortCounts) {
            sortedCountRecord = sortHelper.sort(0,record);
        } else {
            sortedCountRecord = record;
        }
    }

    protected BaseInformationRecords.BaseInformation sortedCountRecord;
    private RecordCountSortHelper sortHelper = new RecordCountSortHelper();

}
