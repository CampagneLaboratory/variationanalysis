package org.campagnelab.dl.genotype.predictions;

import it.unimi.dsi.fastutil.objects.ObjectArraySet;
import it.unimi.dsi.fastutil.objects.ObjectSet;
import org.campagnelab.dl.framework.domains.prediction.Prediction;
import org.campagnelab.dl.genotype.learning.domains.predictions.HomozygousPrediction;
import org.campagnelab.dl.genotype.learning.domains.predictions.SingleGenotypePrediction;

import java.util.Collections;
import java.util.Set;

/**
 * Abstract class to specify genotype predictions
 * Created by rct66 on 12/18/16.
 */
public abstract class AbstractGenotypePrediction extends Prediction {

    /**
     * Genotype called by the model.
     */
    public String predictedGenotype;

    /*
    returns whether or not the prediction is correct, will be used for statistics.
     */
    public abstract boolean isCorrect();

    public Set<String> alleles() {
        return alleles(predictedGenotype);
    }

    public static Set<String> alleles(String genotype) {
        ObjectSet<String> result = new ObjectArraySet<>();
        Collections.addAll(result, genotype.split("[|/]"));
        result.remove("|");
        result.remove("/");
        result.remove("?");
        result.remove(".");
        result.remove("");
        return result;
    }

}
