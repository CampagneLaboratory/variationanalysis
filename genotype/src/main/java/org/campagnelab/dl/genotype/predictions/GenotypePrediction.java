package org.campagnelab.dl.genotype.predictions;

import it.unimi.dsi.fastutil.objects.ObjectArrayList;
import it.unimi.dsi.fastutil.objects.ObjectArraySet;
import it.unimi.dsi.fastutil.objects.ObjectSet;
import org.campagnelab.dl.framework.domains.prediction.Prediction;
import org.campagnelab.dl.genotype.learning.domains.NumDistinctAlleles;
import org.campagnelab.dl.genotype.learning.domains.NumDistinctAllelesInterpreter;
import org.campagnelab.dl.genotype.learning.domains.predictions.HomozygousPrediction;
import org.campagnelab.dl.genotype.learning.domains.predictions.SingleGenotypePrediction;

import java.util.Collections;
import java.util.List;
import java.util.Set;

/**
 * Describes a genotype prediction. Helper method set aggregates individual model predictions.
 * Created by fac2003 on 12/18/16.
 */
public class GenotypePrediction extends AbstractGenotypePrediction{
    private  int predictionIndex;
    /**
     * Genotype called by the model.
     */
    public String calledGenotype;
    public double overallProbability;
    /**
     * True genotype (if available). Useful for performance evaluation.
     */
    public String trueGenotype;
    /**
     * Probability that a genotype is present. Genotypes are indexed using the goby conventions.
     */
    public double probabilityGenotypeCalled[];
    /**
     * Probability that a genotype is not present. Genotypes are indexed using the goby conventions.
     */
    private double probabilityGenotypeNotCalled[];
  //  private HomozygousPrediction homozygousPrediction;
    private SingleGenotypePrediction[] singleGenotypePredictions;
    private int numAlleles;
    private boolean variant;
    public GenotypePrediction() {

    }
    public GenotypePrediction(List<Prediction> predictionList) {
        final Prediction firstPrediction = predictionList.get(0);
        List<Prediction> genoPredList = predictionList.subList(1, 11);

        if (firstPrediction instanceof NumDistinctAlleles) {
            NumDistinctAlleles nda = (NumDistinctAlleles) firstPrediction;
            predictionIndex = nda.index;
            variant = nda.isVariant;
            set(nda, genoPredList.toArray(new SingleGenotypePrediction[genoPredList.size()]));
        }else {
            if (firstPrediction instanceof HomozygousPrediction) {
                HomozygousPrediction homoPred = (HomozygousPrediction) firstPrediction;
                predictionIndex = homoPred.index;
                variant = homoPred.isVariant;
                set(homoPred, genoPredList.toArray(new SingleGenotypePrediction[genoPredList.size()]));
            }
        }
    }


    public void set(HomozygousPrediction homozygousPrediction, SingleGenotypePrediction[] singleGenotypePredictions) {
     //   this.homozygousPrediction = homozygousPrediction;
        this.singleGenotypePredictions = singleGenotypePredictions;
        this.trueGenotype = homozygousPrediction.trueGenotypeFormat;
        if (homozygousPrediction.isHomozygous) {
            calledGenotype = homozygousPrediction.predictedHomozygousGenotype + "/" + homozygousPrediction.predictedHomozygousGenotype;
            overallProbability = homozygousPrediction.probability;
        } else {
            int genotypeIndex = 0;
            double predProbability = 0;
            StringBuffer hetGenotype = new StringBuffer();
            probabilityGenotypeCalled = new double[singleGenotypePredictions.length];
            probabilityGenotypeNotCalled = new double[singleGenotypePredictions.length];
            for (SingleGenotypePrediction singleGenotypePrediction : singleGenotypePredictions) {
                probabilityGenotypeCalled[genotypeIndex] = singleGenotypePrediction.probabilityIsCalled;
                probabilityGenotypeNotCalled[genotypeIndex] = 1 - singleGenotypePrediction.probabilityIsCalled;

                if (singleGenotypePrediction.probabilityIsCalled >= 0.5) {
                    predProbability += singleGenotypePrediction.probabilityIsCalled;
                    numAlleles++;
                    if (hetGenotype.length() > 0) {
                        hetGenotype.append("/");
                    }
                    hetGenotype.append(singleGenotypePrediction.predictedSingleGenotype);
                }
            }
            overallProbability = predProbability / (double) numAlleles;
            calledGenotype = hetGenotype.toString();
        }
    }


    public void set(NumDistinctAlleles numDistinctAlleles, SingleGenotypePrediction[] singleGenotypePredictions) {

        this.singleGenotypePredictions = singleGenotypePredictions;
        this.trueGenotype = numDistinctAlleles.trueGenotypeFormat;
        numAlleles=numDistinctAlleles.predictedValue;

        double predProbability = 0;
        StringBuffer hetGenotype = new StringBuffer();
        ObjectArrayList<SingleGenotypePrediction> list = ObjectArrayList.wrap(singleGenotypePredictions);

        Collections.sort(list, (g1, g2) -> (int) Math.round(1000 * (g2.probabilityIsCalled - g1.probabilityIsCalled)));
        for (SingleGenotypePrediction element : list.subList(0, numAlleles)) {
            if (hetGenotype.length() > 0) {
                hetGenotype.append("/");
            }
            hetGenotype.append(element.predictedSingleGenotype);
            predProbability += element.probabilityIsCalled;
        }
        overallProbability = predProbability / (double) numAlleles;
        calledGenotype = hetGenotype.toString();
    }

    public static Set<String> alleles(String genotype) {
        ObjectSet<String> result = new ObjectArraySet<>();
        Collections.addAll(result, genotype.split("[|/]"));
        return result;
    }

    public Set<String> alleles() {
        return alleles(calledGenotype);
    }

    public boolean isCorrect() {
        Set<String> predictedAlleles = new ObjectArraySet<>();
        for (String s : calledGenotype.split("/")) {
            predictedAlleles.add(s);
        }
        Set<String> trueAlleles = new ObjectArraySet<>();
        for (String s : trueGenotype.split("/")) {
            trueAlleles.add(s);
        }
        Set<String> toIgnore = new ObjectArraySet<String>(new String[]{"?", ".", ""});
        predictedAlleles.removeAll(toIgnore);
        trueAlleles.removeAll(toIgnore);
        return predictedAlleles.equals(trueAlleles);
    }

    public int getPredictionIndex() {
        return predictionIndex;
    }

    public boolean isVariant() {
        return variant;
    }
}
