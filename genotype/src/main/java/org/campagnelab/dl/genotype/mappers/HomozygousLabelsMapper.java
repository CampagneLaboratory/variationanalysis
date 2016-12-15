package org.campagnelab.dl.genotype.mappers;

import org.campagnelab.dl.framework.mappers.MappedDimensions;
import org.campagnelab.dl.somatic.mappers.NoMasksLabelMapper;
import org.campagnelab.dl.varanalysis.protobuf.BaseInformationRecords;
import org.nd4j.linalg.api.ndarray.INDArray;

/**
 * A label that informs if the site is homozygote (label index 11), or indicates the genotype called
 * Created by rct66 on 12/6/16.
 */
public class HomozygousLabelsMapper extends NoMasksLabelMapper<BaseInformationRecords.BaseInformation> {
    @Override
    public int numberOfLabels() {
        return 11;
    }

    int[] indices = new int[]{0, 0};

    @Override
    public void mapLabels(BaseInformationRecords.BaseInformation record, INDArray labels, int indexOfRecord) {
        indices[0] = indexOfRecord;

        for (int labelIndex = 0; labelIndex < 11; labelIndex++) {
            indices[1] = labelIndex;
            labels.putScalar(indices, produceLabel(sortedCountRecord, labelIndex));
        }
    }

    @Override
    public void prepareToNormalize(BaseInformationRecords.BaseInformation record, int indexOfRecord) {
        sortedCountRecord = sortHelper.sort(record);
    }

    private BaseInformationRecords.BaseInformation sortedCountRecord;
    private RecordCountSortHelper sortHelper = new RecordCountSortHelper();

    @Override
    public float produceLabel(BaseInformationRecords.BaseInformation record, int labelIndex) {

        record = sortedCountRecord;
        if (!getHomozygous(record)) {
            // this site is heterozygous. The Allele will only be encoded in the other outputs.
            return (labelIndex == 10) ? 1 : 0;
        } else {
            // the site is homozygous. The allele is encoded here:
            if (labelIndex >= record.getSamples(0).getCountsCount()) {
                return 0;
            } else {
                return record.getSamples(0).getCounts(labelIndex).getIsCalled() ? 1 : 0;
            }
        }
    }


    private boolean getHomozygous(BaseInformationRecords.BaseInformation record) {
        record = sortedCountRecord;
        //count number of called alleles
        int numAlleles = 0;
        for (int i = 0; i < record.getSamples(0).getCountsCount(); i++) {
            if (record.getSamples(0).getCounts(i).getIsCalled()) {
                numAlleles++;
            }
        }
        return (numAlleles == 1);
    }


    @Override
    public MappedDimensions dimensions() {
        return new MappedDimensions(numberOfLabels());
    }


}
