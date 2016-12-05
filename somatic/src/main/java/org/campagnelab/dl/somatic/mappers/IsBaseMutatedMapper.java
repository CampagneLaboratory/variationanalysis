package org.campagnelab.dl.somatic.mappers;

import org.campagnelab.dl.varanalysis.protobuf.BaseInformationRecords;
import org.nd4j.linalg.api.ndarray.INDArray;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;

/**
 * A label with 6 floats: first is probability that site is not mutated. Next floats are probability
 * that a genotype (sorted by decreasing count) is a somatic mutation.
 * Created by fac2003 on 5/12/2016.
 */
public class IsBaseMutatedMapper extends NoMasksLabelMapper<BaseInformationRecords.BaseInformation> {
    int[] indices = new int[]{0, 0};

    private float[] labels = new float[numberOfLabels()];

    private ArrayList<BaseInformationRecords.CountInfo> sortedCounts;

    @Override
    public void prepareToNormalize(BaseInformationRecords.BaseInformation record, int indexOfRecord) {
        int numLabels = numberOfLabels();
        final List<BaseInformationRecords.CountInfo> counts = record.getSamples(record.getSamplesCount() - 1).getCountsList();
        ArrayList<BaseInformationRecords.CountInfo> sorted = new ArrayList<>();
        sorted.addAll(counts);
        Collections.sort(sorted, (o1, o2) ->
                (o2.getGenotypeCountForwardStrand() + o2.getGenotypeCountReverseStrand()) - (o1.getGenotypeCountForwardStrand() + o1.getGenotypeCountReverseStrand())
        );
        sortedCounts = sorted;
        Arrays.fill(labels, 0f);
        labels[0] = record.getMutated() ? 0 : 1;
        if (record.getMutated()) {
            for (int i = 1; i < numLabels; i++) {
                labels[i] = sorted.get(i - 1).getToSequence().equals(record.getMutatedBase()) ? 1 : 0;
            }
        }
    /* labels must sum to 1.
    float sum = 0;
        for (int i = 0; i < numLabels; i++) {
            sum += labels[i];
        }
        System.out.println("sum labels=" + sum);
   */
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
    public int numberOfLabels() {
        return AbstractFeatureMapper.MAX_GENOTYPES + 1;
    }


    @Override
    public float produceLabel(BaseInformationRecords.BaseInformation record, int labelIndex) {
        return labels[labelIndex];
    }

}
