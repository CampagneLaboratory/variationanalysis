package org.campagnelab.dl.genotype.predictions;

import it.unimi.dsi.fastutil.objects.ObjectArrayList;
import org.campagnelab.dl.framework.domains.prediction.Prediction;
import org.campagnelab.dl.genotype.learning.domains.NumDistinctAllelesOutputLayerPrediction;
import org.campagnelab.dl.genotype.learning.domains.predictions.SingleGenotypePrediction;
import org.campagnelab.dl.varanalysis.protobuf.BaseInformationRecords;

import java.util.Collections;
import java.util.List;

/**
 * A prediction for models that encode the number of distinct alleles.
 */
public class NumDistinctAlleleGenotypePrediction extends GenotypePrediction {

    private final BaseInformationRecords.BaseInformation record;
    private boolean useNumAlleles = false;

    public NumDistinctAlleleGenotypePrediction(BaseInformationRecords.BaseInformation record, double decisionThreshold, List<Prediction> predictionList) {
        this.DECISION_THRESHOLD = decisionThreshold;
        this.record = record;
        set(predictionList);
    }

    public void set(List<Prediction> predictionList) {
        final Prediction firstPrediction = predictionList.get(0);
        List<Prediction> genoPredList = predictionList.subList(1, 11);
        NumDistinctAllelesOutputLayerPrediction nda = (NumDistinctAllelesOutputLayerPrediction) firstPrediction;
        index = nda.index;
        set(nda, genoPredList.toArray(new SingleGenotypePrediction[genoPredList.size()]),
                (MetadataPrediction) predictionList.get(11));

    }

    private double DECISION_THRESHOLD = 0.5d;

    public void set(NumDistinctAllelesOutputLayerPrediction numDistinctAlleles, SingleGenotypePrediction[] singleGenotypePredictions, MetadataPrediction metaData) {
        double threshold = DECISION_THRESHOLD;
        ObjectArrayList<SingleGenotypePrediction> list = ObjectArrayList.wrap(singleGenotypePredictions);
        list.sort((g1, g2) -> Double.compare(g2.probabilityIsCalled, g1.probabilityIsCalled));

        int numAlleles = useNumAlleles ? numDistinctAlleles.predictedValue : calculateNumAlleles(threshold, list);

        double predProbability = numDistinctAlleles.probability;

        StringBuffer hetGenotype = new StringBuffer();
        isPredictedIndel = false;
        int indexFirstNotCalled = numAlleles + 1;
        // will hold the probability of the called allele with lowest probability
        double weakestProbabilityCalled = Double.POSITIVE_INFINITY;
        for (SingleGenotypePrediction element : list.subList(0, numAlleles)) {
            if (hetGenotype.length() > 0) {
                hetGenotype.append("/");
            }
            hetGenotype.append(element.predictedSingleGenotype(record, metaData));
            final boolean alleleIsIndel = element.isPredictedIndel(record, metaData);
            this.isPredictedIndel |= alleleIsIndel;
            predProbability += Math.max(element.probabilityIsCalled, 1 - element.probabilityIsCalled);
            weakestProbabilityCalled = Math.min(weakestProbabilityCalled, element.probabilityIsCalled);
        }
        this.trueGenotype = extractTrueGenotype(record, singleGenotypePredictions, metaData);
        double largestProbabilityNotCalled = list.get(indexFirstNotCalled).probabilityIsCalled;
        double confidenceInPrediction = weakestProbabilityCalled - largestProbabilityNotCalled;
        overallProbability = confidenceInPrediction;
        this.isVariantProbability = overallProbability;
        predictedGenotype = hetGenotype.toString();
        this.isIndel = metaData.isIndel;
        this.isVariant = metaData.isVariant;
        this.referenceGobyIndex = metaData.referenceGobyIndex;
    }

    private int calculateNumAlleles(double threshold, ObjectArrayList<SingleGenotypePrediction> list) {
        int numAlelles = 0;
        for (SingleGenotypePrediction genotypePrediction : list) {
            if (genotypePrediction.probabilityIsCalled >= threshold) {
                numAlelles += 1;
            } else {
                break;
            }
        }
        return numAlelles;
    }

    private String extractTrueGenotype(BaseInformationRecords.BaseInformation record, SingleGenotypePrediction[] singleGenotypePredictions, MetadataPrediction metaData) {
        String result = "";
        for (SingleGenotypePrediction single : singleGenotypePredictions) {

            if (single.trueIsCalled) {
                if (result.length() > 0 && !result.endsWith("/")) {
                    result += "/";
                }
                result += single.predictedSingleGenotype(record, metaData);
            }
        }
        return result;
    }


}