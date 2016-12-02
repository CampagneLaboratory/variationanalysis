package org.campagnelab.dl.genotype.mappers;

import org.campagnelab.dl.somatic.mappers.NoMasksLabelMapper;
import org.campagnelab.dl.varanalysis.protobuf.BaseInformationRecords;
import org.nd4j.linalg.api.ndarray.INDArray;

/**
 * Label: frequency of somatic mutation.
 * Created by fac2003 on 11/8/16.
 */
public class GenotypeLabelMapper extends NoMasksLabelMapper<BaseInformationRecords.BaseInformation> {
    @Override
    public int numberOfLabels() {
        return 2;
    }

    int[] indices = new int[]{0, 0};
    int genotypeIndex;

    public GenotypeLabelMapper(int genotypeIndex){
        this.genotypeIndex = genotypeIndex;
    }

    @Override
    public void mapLabels(BaseInformationRecords.BaseInformation record, INDArray labels, int indexOfRecord) {
        indices[0] = indexOfRecord;
        for (int labelIndex = 0; labelIndex < numberOfLabels(); labelIndex++) {
            indices[1] = labelIndex;
            labels.putScalar(indices, produceLabel(record, labelIndex));
        }
    }

    @Override
    public float produceLabel(BaseInformationRecords.BaseInformation record, int labelIndex) {
        assert labelIndex == 0 || labelIndex == 1 : "only one label.";
        boolean isCalled;
        if (genotypeIndex >= record.getSamples(0).getCountsCount()){
            isCalled = false;
        } else {
            isCalled = record.getSamples(0).getCounts(genotypeIndex).getIsCalled();
        }
        if (labelIndex == 0) {
            // first index is 1 when site is  mutated.
            return isCalled?1:0;
        } else {
            // second index is 1 when site is not mutated.

            return !isCalled?1:0;
        }

    }




}
