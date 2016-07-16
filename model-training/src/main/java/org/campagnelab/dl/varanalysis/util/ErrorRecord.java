package org.campagnelab.dl.varanalysis.util;

import org.deeplearning4j.nn.multilayer.MultiLayerNetwork;
import org.nd4j.linalg.api.ndarray.INDArray;

import java.util.Comparator;

/**
 * Created by fac2003 on 7/15/16.
 */
public class ErrorRecord {
    /* A measure of how wrong this prediction was. Larger is more wronger.
         */
    public float wrongness;
    public INDArray features;

    public INDArray label;

    public ErrorRecord(float wrongness, INDArray features, INDArray label) {
        this.wrongness = wrongness;
        this.features = features;
        this.label = label;
    }

    @Override
    public int hashCode() {
        return this.features.hashCode() ^ this.label.hashCode();
    }

    @Override
    public boolean equals(Object obj) {
        if (!(obj instanceof ErrorRecord)) {
            return false;
        }
        ErrorRecord errorRecord = (ErrorRecord) obj;
        return features.equals(errorRecord.features) && label.equals(errorRecord.label);
    }

    public static Comparator<ErrorRecord> INCREASING_SCORE_COMPARATOR = (a, b) -> {

        if (a.equals(b)) {
            return 0;
        }
        return Float.compare(a.wrongness, b.wrongness);
    };
    public static Comparator<ErrorRecord> DECREASING_SCORE_COMPARATOR = (a, b) -> {

        if (a.equals(b)) {
            return 0;
        }
        return Float.compare( b.wrongness,a.wrongness);
    };

    public void updateWrongness(INDArray features, MultiLayerNetwork net) {

        INDArray predictedLabels = net.output(features, false);
        this.wrongness = ErrorRecord.calculateWrongness(0, predictedLabels, label);
    }

    @Override
    public String toString() {
        return String.format("%f features=%s labels=%s", wrongness, features, label);
    }

    public static boolean isWrongPrediction(int exampleIndex, INDArray predictedLabels, INDArray labels) {

        final int positiveLabelMutated = 0;
        final int negativeLabel = 1;
        final boolean predictedPositive = predictedLabels.getDouble(exampleIndex, positiveLabelMutated) > predictedLabels.getDouble(exampleIndex, negativeLabel);
        final double trueLabelPositive = labels.getDouble(exampleIndex, positiveLabelMutated);

        return predictedPositive && trueLabelPositive == 0 ||
                !predictedPositive && trueLabelPositive > 0;

    }

    public static float calculateWrongness(int exampleIndex, INDArray predictedLabels, INDArray labels) {
        final int positiveLabelMutated = 0;
        final int negativeLabel = 1;
        final boolean predictedPositive = predictedLabels.getDouble(exampleIndex, positiveLabelMutated) > predictedLabels.getDouble(exampleIndex, negativeLabel);
        if (predictedPositive) {
            return predictedLabels.getFloat(exampleIndex, positiveLabelMutated);
        } else {
            return predictedLabels.getFloat(exampleIndex, negativeLabel);
        }
    }
}
