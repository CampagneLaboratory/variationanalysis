package org.campagnelab.dl.genotype.learning.domains.predictions;

import org.campagnelab.dl.framework.domains.prediction.Prediction;
import org.campagnelab.dl.genotype.predictions.MetadataPrediction;
import org.campagnelab.dl.varanalysis.protobuf.BaseInformationRecords;

/**
 * Stores information about the call for a single allele.
 * Created by rct66 on 11/12/16.
 */
public class SingleGenotypePrediction extends Prediction {

    /**
     * The allele called/predicted to be present in the sample. e.g., "A" or "A--"
     */
    public String predictedSingleGenotype(BaseInformationRecords.BaseInformationOrBuilder record, MetadataPrediction metaData) {
        if (sortedCountIndex>= metaData.sorted2OriginalCountIndices.length) {
            // ONLY 3 sorted genotypes.
            return ".";
        }
        final int originalCountIndex = metaData.sorted2OriginalCountIndices[sortedCountIndex];
        if (record == null) {
            // no record in training phase. We return an integer with the index of the original goby count
            // (e.g., "5" for the first indel).
            return Integer.toString(originalCountIndex);
        } else {
            final BaseInformationRecords.SampleInfo sample = record.getSamples(0);
            if (originalCountIndex >= sample.getCountsCount()) {
                return ".";
            } else {
                return sample.getCountsOrBuilder(originalCountIndex).getToSequence();
            }
        }
    }

    /**
     * The probability the allele is present, given the model and the data.
     */
    public double probabilityIsCalled;
    /**
     * Whether the allele is present in the ground-truth data.
     */
    public boolean trueIsCalled;

    /**
     * Index of the sorted count corresponding to this genotytpe.
     */
    public int sortedCountIndex;

    public boolean isPredictedIndel(BaseInformationRecords.BaseInformation record, MetadataPrediction metaData) {
        return sortedCountIndex < 3 && (metaData.sorted2OriginalCountIndices[sortedCountIndex] > 4);
    }
}

