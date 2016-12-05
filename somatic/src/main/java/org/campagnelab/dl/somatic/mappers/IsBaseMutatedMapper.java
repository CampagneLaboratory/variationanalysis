package org.campagnelab.dl.somatic.mappers;

import it.unimi.dsi.fastutil.objects.ObjectArrayList;
import it.unimi.dsi.fastutil.objects.ObjectSet;
import org.campagnelab.dl.framework.iterators.ConcatFeatureMapper;
import org.campagnelab.dl.framework.iterators.ConcatLabelMapper;
import org.campagnelab.dl.framework.mappers.OneHotBaseFeatureMapper;
import org.campagnelab.dl.varanalysis.protobuf.BaseInformationRecords;
import org.nd4j.linalg.api.ndarray.INDArray;

import java.util.Arrays;
import java.util.function.Function;

/**
 * Created by fac2003 on 11/8/16.
 */
public class IsBaseMutatedMapper extends NoMasksLabelMapper<BaseInformationRecords.BaseInformation> {
    private final ConcatLabelMapper<BaseInformationRecords.BaseInformation> delegate=null;
    int[] indices = new int[]{0, 0};

    ObjectArrayList<String> genotypes;

    @Override
    public void prepareToNormalize(BaseInformationRecords.BaseInformation record, int indexOfRecord) {


    }
    public IsBaseMutatedMapper(int contextSize, Function<BaseInformationRecords.BaseInformationOrBuilder, String> function) {
   //     GenomicContextMapper<BaseInformationRecords.BaseInformationOrBuilder>[] refContext = new OneHotBaseFeatureMapper[contextSize];
     //   for (int i = 0; i < contextSize; i++) {
       //     refContext[i] = new OneHotBaseFeatureMapper<>(i, function);
        //}
        //delegate = new ConcatLabelMapper<>(refContext);
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
        assert labelIndex == 0 || labelIndex == 1 : "only one label.";

        // first index is 1 when site is  mutated.
        if (labelIndex == 0) return record.getMutated() ?1 : 0;
        // second index is 1 when site is not mutated.
        return record.getMutated() ? 0 : 1;
    }

    public int indexOf(String genotype) {
return 0;
    }
}
