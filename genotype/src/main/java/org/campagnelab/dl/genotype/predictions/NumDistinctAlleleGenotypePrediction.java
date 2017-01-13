package org.campagnelab.dl.genotype.predictions;

import it.unimi.dsi.fastutil.objects.ObjectArrayList;
import org.campagnelab.dl.framework.domains.prediction.Prediction;
import org.campagnelab.dl.genotype.learning.domains.NumDistinctAllelesOutputLayerPrediction;
import org.campagnelab.dl.genotype.learning.domains.predictions.SingleGenotypePrediction;

import java.util.Collections;
import java.util.List;

/**
 * A prediction for models that encode the number of distinct alleles.
 */
public class NumDistinctAlleleGenotypePrediction extends GenotypePrediction {

    private boolean useNumAlleles = false;

    public NumDistinctAlleleGenotypePrediction(List<Prediction> predictionList) {
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

    public void set(NumDistinctAllelesOutputLayerPrediction numDistinctAlleles, SingleGenotypePrediction[] singleGenotypePredictions, MetadataPrediction metaData) {
        double threshold = 0.5;
        ObjectArrayList<SingleGenotypePrediction> list = ObjectArrayList.wrap(singleGenotypePredictions);
        list.sort((g1, g2) -> Double.compare(g2.probabilityIsCalled, g1.probabilityIsCalled));

        int numAlleles = useNumAlleles ? numDistinctAlleles.predictedValue : calculateNumAlleles(threshold, list);

        double predProbability = numDistinctAlleles.probability;

        StringBuffer hetGenotype = new StringBuffer();

        for (SingleGenotypePrediction element : list.subList(0, numAlleles)) {
            if (hetGenotype.length() > 0) {
                hetGenotype.append("/");
            }
            hetGenotype.append(element.predictedSingleGenotype);
            predProbability += Math.max(element.probabilityIsCalled, 1 - element.probabilityIsCalled);
        }
        this.trueGenotype = extractTrueGenotype(singleGenotypePredictions);
        overallProbability = predProbability / (double) (numAlleles + 1);
        this.isVariantProbability = overallProbability;
        predictedGenotype = hetGenotype.toString();
        this.isIndel = metaData.isIndel;
        this.isVariant = metaData.isVariant;
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

    private String extractTrueGenotype(SingleGenotypePrediction[] singleGenotypePredictions) {
        String result = "";
        for (SingleGenotypePrediction single : singleGenotypePredictions) {

            if (single.trueIsCalled) {
                if (result.length() > 0 && !result.endsWith("/")) {
                    result += "/";
                }
                result += single.predictedSingleGenotype;
            }
        }
        return result;
    }


}