package org.campagnelab.dl.genotype.mappers;

import org.campagnelab.dl.framework.mappers.MappedDimensions;
import org.campagnelab.dl.somatic.mappers.NoMasksLabelMapper;
import org.campagnelab.dl.varanalysis.protobuf.BaseInformationRecords;
import org.nd4j.linalg.api.ndarray.INDArray;

/**
 * Label: frequency of somatic mutation.
 * Created by rct66 on 12/6/16.
 */
public class GenotypeLabelsMapper extends NoMasksLabelMapper<BaseInformationRecords.BaseInformation> {
    @Override
    public int numberOfLabels() {
        return 2;
    }

    int genotypeIndex;
    int[] indices = new int[]{0, 0};



    public GenotypeLabelsMapper(int genotypeIndex){
        this.genotypeIndex = genotypeIndex;
    }

    @Override
    public void mapLabels(BaseInformationRecords.BaseInformation record, INDArray labels, int indexOfRecord) {
        indices[0] = indexOfRecord;
        for (int labelIndex = 0; labelIndex < 2; labelIndex++) {
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


    @Override
    public MappedDimensions dimensions() {
        return new MappedDimensions(new int[]{2});
    }


}
