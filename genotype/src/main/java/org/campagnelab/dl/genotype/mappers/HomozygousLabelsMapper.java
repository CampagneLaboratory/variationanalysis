package org.campagnelab.dl.genotype.mappers;

import org.campagnelab.dl.framework.mappers.MappedDimensions;
import org.campagnelab.dl.somatic.mappers.NoMasksLabelMapper;
import org.campagnelab.dl.varanalysis.protobuf.BaseInformationRecords;
import org.nd4j.linalg.api.ndarray.INDArray;

/**
 * Label: frequency of somatic mutation.
 * Created by rct66 on 12/6/16.
 */
public class HomozygousLabelsMapper extends NoMasksLabelMapper<BaseInformationRecords.BaseInformation> {
    @Override
    public int numberOfLabels() {
        return 11;
    }

    int[] indices = new int[]{0, 0};


    private Boolean isHomozygous;

    @Override
    public void mapLabels(BaseInformationRecords.BaseInformation record, INDArray labels, int indexOfRecord) {
        indices[0] = indexOfRecord;
        for (int labelIndex = 0; labelIndex < 11; labelIndex++) {
            indices[1] = labelIndex;
            labels.putScalar(indices, produceLabel(record, labelIndex));
        }
    }

    @Override
    public float produceLabel(BaseInformationRecords.BaseInformation record, int labelIndex) {

        if (!getHomozygous(record)) {
            return (labelIndex == 10) ? 1 : 0;
        }
        if (labelIndex >= record.getSamples(0).getCountsCount()){
            return 0;
        } else {
            return record.getSamples(0).getCounts(labelIndex).getIsCalled()?1:0;
        }
    }


    private boolean getHomozygous(BaseInformationRecords.BaseInformation record) {
        if (isHomozygous != null) {
            return isHomozygous;
        }
        //count number of called alleles
        int numAlleles = 0;
        for (int i = 0; i < record.getSamples(0).getCountsCount(); i++){
            if (record.getSamples(0).getCounts(i).getIsCalled()){
                numAlleles++;
            }
        }
        isHomozygous = (numAlleles == 1);
        return isHomozygous;
    }



    @Override
    public MappedDimensions dimensions() {
        return new MappedDimensions(new int[]{11});
    }


}
