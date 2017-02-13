package org.campagnelab.dl.genotype.mappers;

import org.campagnelab.dl.framework.mappers.FeatureNameMapper;
import org.campagnelab.dl.framework.mappers.NoMaskFeatureMapper;
import org.campagnelab.dl.varanalysis.protobuf.BaseInformationRecords;
import org.nd4j.linalg.api.ndarray.INDArray;

import java.util.Arrays;

/**
 * Determine the number of genotypes with at least 90% base support, which do not match the reference.
 * Created by fac2003 on 1/14/17.
 */
public class NaiveNumAlleleMapper<T>
        extends NoMaskFeatureMapper<BaseInformationRecords.BaseInformationOrBuilder>
        implements FeatureNameMapper<BaseInformationRecords.BaseInformationOrBuilder> {
    @Override
    public int numberOfFeatures() {
        return MAX_GENOTYPES;
    }

    public NaiveNumAlleleMapper(int sampleIndex) {
        this.sampleIndex = sampleIndex;
    }

    int sampleIndex = 0;
    int counts[] = new int[MAX_GENOTYPES];
    private static int MAX_GENOTYPES = 3;
    int numAlleles = 0;

    @Override
    public void prepareToNormalize(BaseInformationRecords.BaseInformationOrBuilder record, int indexOfRecord) {
        int sumCounts = 0;
        Arrays.fill(counts, 0);
        for (int i = 0; i < MAX_GENOTYPES; i++) {
            BaseInformationRecords.CountInfo counts = record.getSamples(sampleIndex).getCounts(i);
            if (!counts.getMatchesReference()) {
                this.counts[i] += counts.getGenotypeCountForwardStrand() + counts.getGenotypeCountReverseStrand();
                sumCounts += this.counts[i];
            }
        }
        numAlleles = 0;
        int cumulativeCount = 0;
        for (int i = 0; i < MAX_GENOTYPES; i++) {
            cumulativeCount += this.counts[i];
            if (cumulativeCount >= sumCounts * 0.9) {
                numAlleles = i + 1;
                break;
            }
        }
     //   System.out.println("num Alleles:" + numAlleles);
    }

    private static final int[] indices = new int[]{0, 0};

    @Override
    public void mapFeatures(BaseInformationRecords.BaseInformationOrBuilder record, INDArray inputs, int indexOfRecord) {

        indices[0] = indexOfRecord;
        for (int featureIndex = 0; featureIndex < numberOfFeatures(); featureIndex++) {
            indices[1] = featureIndex;
            inputs.putScalar(indices, produceFeature(record, featureIndex));
        }
    }

    @Override
    public float produceFeature(BaseInformationRecords.BaseInformationOrBuilder record, int featureIndex) {
        return featureIndex == numAlleles ? 1 : 0;
    }

    @Override
    public String getFeatureName(int featureIndex) {
        return "NaiveNumAlleleMapper" + featureIndex;
    }
}
