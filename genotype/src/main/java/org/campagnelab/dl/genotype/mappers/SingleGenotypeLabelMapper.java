package org.campagnelab.dl.genotype.mappers;

import org.campagnelab.dl.framework.mappers.MappedDimensions;
import org.campagnelab.dl.somatic.mappers.NoMasksLabelMapper;
import org.campagnelab.dl.varanalysis.protobuf.BaseInformationRecords;
import org.nd4j.linalg.api.ndarray.INDArray;

/**
 * Label: whether a genotype index (original Goby index, before sorting) is called or not.
 * Created by rct66 on 12/6/16.
 */
public class SingleGenotypeLabelMapper extends NoMasksLabelMapper<BaseInformationRecords.BaseInformation> {
    private final boolean sortCounts;
    private final float epsilon;
    private int[] indexPermutation;

    @Override
    public int numberOfLabels() {
        return 2;
    }

    int genotypeIndex;
    int[] indices = new int[]{0, 0};

    public SingleGenotypeLabelMapper(int genotypeIndex, boolean sort) {
        this(genotypeIndex, sort, 0);
    }

    public SingleGenotypeLabelMapper(int genotypeIndex, boolean sort, float epsilon) {
        this.genotypeIndex = genotypeIndex;
        this.sortCounts = sort;
        this.epsilon = epsilon;
    }

    @Override
    public void mapLabels(BaseInformationRecords.BaseInformation record, INDArray labels, int indexOfRecord) {
        indices[0] = indexOfRecord;
        for (int labelIndex = 0; labelIndex < 2; labelIndex++) {
            indices[1] = labelIndex;
            labels.putScalar(indices, produceLabel(sortedCountRecord, labelIndex));
        }
    }

    private BaseInformationRecords.BaseInformation sortedCountRecord;
    private RecordCountSortHelper sortHelper = new RecordCountSortHelper();

    @Override
    public float produceLabel(BaseInformationRecords.BaseInformation record, int labelIndex) {
        assert labelIndex == 0 || labelIndex == 1 : "only one label.";
        boolean isCalled;
        record = sortedCountRecord;
        if (genotypeIndex >= record.getSamples(0).getCountsCount()) {
            isCalled = false;
        } else {
            isCalled = record.getSamples(0).getCounts(indexPermutation[genotypeIndex]).getIsCalled();
        }
        if (labelIndex == 0) {
            // first index is 1 when site is  called.
            return isCalled ? 1 - epsilon : epsilon;
        } else {
            // second index is 1 when site is not called.
            return !isCalled ? 1 - epsilon : epsilon;
        }
    }

    @Override
    public void prepareToNormalize(BaseInformationRecords.BaseInformation record, int indexOfRecord) {
        int sortedIndex = 0;
        indexPermutation = new int[sortedCountRecord.getSamples(0).getCountsCount()];
        if (sortCounts) {
            sortedCountRecord = sortHelper.sort(record);

            for (BaseInformationRecords.CountInfo count : sortedCountRecord.getSamples(0).getCountsList()) {
                indexPermutation[count.getGobyGenotypeIndex()] = sortedIndex++;
            }
        } else {
            sortedCountRecord = record;
            for (BaseInformationRecords.CountInfo count : sortedCountRecord.getSamples(0).getCountsList()) {
                indexPermutation[sortedIndex] = sortedIndex;
                sortedIndex++;
            }
        }

    }


    @Override
    public MappedDimensions dimensions() {
        return new MappedDimensions(2);
    }


}
