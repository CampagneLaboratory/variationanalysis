package org.campagnelab.dl.genotype.learning.domains;

import org.campagnelab.dl.somatic.mappers.NoMasksLabelMapper;
import org.campagnelab.dl.varanalysis.protobuf.BaseInformationRecords;
import org.nd4j.linalg.api.ndarray.INDArray;

import java.util.function.Predicate;

/**
 * A functional label mapper for boolean values.
 * Created by fac2003 on 12/23/16.
 */
public class BooleanLabelMapper<RecordType> extends NoMasksLabelMapper<RecordType> {
    private Predicate<RecordType> predicate;

    public BooleanLabelMapper(Predicate<RecordType> predicate) {
        this.predicate = predicate;
    }

    @Override
    public int numberOfLabels() {
        return 2;
    }

    int[] indices = new int[]{0, 0};

    @Override
    public void mapLabels(RecordType record, INDArray labels, int indexOfRecord) {
        indices[0] = indexOfRecord;

        for (int labelIndex = 0; labelIndex < numberOfLabels(); labelIndex++) {
            indices[1] = labelIndex;
            labels.putScalar(indices, produceLabel(record, labelIndex));
        }
    }


    @Override
    public float produceLabel(RecordType record, int labelIndex) {
        switch (labelIndex) {
            case 0:
                return predicate.test(record) ? 1 : 0;
            case 1:
                return !predicate.test(record) ? 1 : 0;
            default:
                throw new RuntimeException("Label index out of range, must be 0 or 1: " + labelIndex);
        }

    }
}
