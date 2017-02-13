package org.campagnelab.dl.genotype.mappers;

import org.campagnelab.dl.framework.mappers.LabelMapper;
import org.campagnelab.dl.framework.mappers.MappedDimensions;
import org.campagnelab.dl.genotype.helpers.GenotypeHelper;
import org.campagnelab.dl.somatic.mappers.NoMasksLabelMapper;
import org.campagnelab.dl.varanalysis.protobuf.BaseInformationRecords;
import org.nd4j.linalg.api.ndarray.INDArray;

/**
 * Stores meta-data in a virtual label, not used for learning.
 * Created by fac2003 on 12/22/16.
 */
public class MetaDataLabelMapper extends NoMasksLabelMapper<BaseInformationRecords.BaseInformation> {

    public static final int NUM_LABELS = 2;
    public static final int IS_VARIANT_FEATURE_INDEX = 0;
    public static final int IS_INDEL_FEATURE_INDEX = 1;

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
            labels.putScalar(indices, produceLabel(record, labelIndex));
        }
    }


    @Override
    public float produceLabel(BaseInformationRecords.BaseInformation record, int labelIndex) {
        switch (labelIndex) {
            case IS_VARIANT_FEATURE_INDEX:
                return record.getSamples(0).getIsVariant() ? 1 : 0;
            case IS_INDEL_FEATURE_INDEX:
                final String trueGenotype = record.getTrueGenotype();
                return GenotypeHelper.isIndel(record.getReferenceBase(), trueGenotype) ? 1 : 0;
            default:
                throw new RuntimeException("No such labelIndex: " + labelIndex);
        }
    }


    @Override
    public MappedDimensions dimensions() {
        return new MappedDimensions(numberOfLabels());
    }
}

