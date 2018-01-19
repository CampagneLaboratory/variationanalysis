package org.campagnelab.dl.genotype.mappers;

import org.campagnelab.dl.framework.mappers.MappedDimensions;
import org.campagnelab.dl.varanalysis.protobuf.BaseInformationRecords;
import org.nd4j.linalg.api.ndarray.INDArray;

/**
 * Label: whether a genotype index (original Goby index, before sorting) is called or not.
 * Created by rct66 on 12/6/16.
 */
public class SingleGenotypeLabelMapper extends RecordCountSortingLabelMapperImpl {

    @Override
    public int numberOfLabels() {
        return 2;
    }

    private int sortedGenotypeIndex;
    private int[] indices = new int[]{0, 0};
    private final float epsilon;

    public SingleGenotypeLabelMapper(int sortedGenotypeIndex, boolean sort) {
        this(sortedGenotypeIndex, sort, 0);
    }

    public SingleGenotypeLabelMapper(int sortedGenotypeIndex, boolean sort, float epsilon) {
        super(sort);
        this.sortedGenotypeIndex = sortedGenotypeIndex;
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

    @Override
    public float produceLabel(BaseInformationRecords.BaseInformation record, int labelIndex) {
        assert labelIndex == 0 || labelIndex == 1 : "only one label.";
        boolean isCalled;
        record = sortedCountRecord;
        if (sortedGenotypeIndex >= record.getSamples(this.sampleIndex).getCountsCount() || sortedGenotypeIndex>= GenotypeFeatureMapper.MAX_GENOTYPES) {
            isCalled = false;
        } else {
            isCalled = record.getSamples(this.sampleIndex).getCounts(sortedGenotypeIndex).getIsCalled();
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
    public MappedDimensions dimensions() {
        return new MappedDimensions(2);
    }


}
