package org.campagnelab.dl.genotype.mappers;

import org.campagnelab.dl.framework.mappers.MappedDimensions;
import org.campagnelab.dl.somatic.mappers.NoMasksLabelMapper;
import org.campagnelab.dl.varanalysis.protobuf.BaseInformationRecords;
import org.nd4j.linalg.api.ndarray.INDArray;

/**
 * A full softmax mapper that offers homozygous first, heterozygous first/second, and homozygous second no-call.
 * Will eventually be extendible to more alternates and more alleles per genotype (ie non-diploid plants)
 *
 * REQUIRES sorted features
 *
 * Created by rct66 on 12/20/16.
 *
 *
 */
public class CombinedLabelsMapper extends NoMasksLabelMapper<BaseInformationRecords.BaseInformation> {

    public static final int NUM_LABELS = 4;
private float epsilon=0;

    public CombinedLabelsMapper(float epsilon) {
        this.epsilon = epsilon;
    }

    @Override
    public int numberOfLabels() {
        return NUM_LABELS;
    }


    int[] indices = new int[]{0, 0};

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
        sortedCountRecord = sortHelper.sort(record);
    }

    private BaseInformationRecords.BaseInformation sortedCountRecord;
    private RecordCountSortHelper sortHelper = new RecordCountSortHelper();


    @Override
    public float produceLabel(BaseInformationRecords.BaseInformation record, int labelIndex) {
        int correctLabelIndex;
        boolean firstCalled = sortedCountRecord.getSamples(0).getCounts(0).getIsCalled();
        boolean secondCalled = sortedCountRecord.getSamples(0).getCounts(1).getIsCalled();
        boolean otherCalled = false;
        for (int i = 2; i < sortedCountRecord.getSamples(0).getCountsCount(); i++){
            otherCalled |= sortedCountRecord.getSamples(0).getCounts(i).getIsCalled();
        }
        if (otherCalled){
            //other call case
            correctLabelIndex = 3;
        } else if (firstCalled){
            if (secondCalled){
                //both case (het)
                correctLabelIndex = 1;
            } else {
                //homozyg most counts
                correctLabelIndex = 0;
            }
        } else if (secondCalled){
            //homozyg smaller counts
            correctLabelIndex = 2;
        } else {
            //no call case
            correctLabelIndex = 3;
        }
        return (labelIndex == correctLabelIndex)?1-epsilon:epsilon/(numberOfLabels()-1);
    }




    @Override
    public MappedDimensions dimensions() {
        return new MappedDimensions(numberOfLabels());
    }


}
