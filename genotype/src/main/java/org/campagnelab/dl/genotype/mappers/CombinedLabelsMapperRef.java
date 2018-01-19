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
public class CombinedLabelsMapperRef extends NoMasksLabelMapper<BaseInformationRecords.BaseInformation> {

    public static final int NUM_LABELS = 4;
    private float epsilon=0;

    public CombinedLabelsMapperRef(float epsilon) {
        this.epsilon = epsilon;
    }

    public CombinedLabelsMapperRef() {
        this(0);
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
        boolean firstCalled = sortedCountRecord.getSamples(this.sampleIndex).getCounts(0).getIsCalled();
        boolean secondCalled = sortedCountRecord.getSamples(this.sampleIndex).getCounts(1).getIsCalled();
        boolean firstIsRef = sortedCountRecord.getSamples(this.sampleIndex).getCounts(0).getMatchesReference();
        boolean secondIsRef = sortedCountRecord.getSamples(this.sampleIndex).getCounts(1).getMatchesReference();;
        boolean otherCalled = false;

        for (int i = 2; i < sortedCountRecord.getSamples(this.sampleIndex).getCountsCount(); i++){
            otherCalled |= sortedCountRecord.getSamples(this.sampleIndex).getCounts(i).getIsCalled();
        }
        if (otherCalled){
            //other call case
            correctLabelIndex = 3;
        } else if (firstCalled && secondCalled) {
            //both called
            if (firstIsRef || secondIsRef){
                correctLabelIndex = 1;
            } else {
                correctLabelIndex = 3;
            }
        } else if (firstCalled){
            //first is called but no second
            if (firstIsRef){
                correctLabelIndex = 0;
            } else {
                correctLabelIndex = 2;
            }
        } else if (secondCalled){
            //second is called but not first
            if (secondIsRef){
                correctLabelIndex = 0;
            } else {
                correctLabelIndex = 2;
            }
        } else {
            //no calls at all
            correctLabelIndex = 3;
        }

        return (labelIndex == correctLabelIndex)?1-epsilon:epsilon/(numberOfLabels()-1);
    }




    @Override
    public MappedDimensions dimensions() {
        return new MappedDimensions(numberOfLabels());
    }


}
