package org.campagnelab.dl.genotype.learning.domains.predictions;

import org.campagnelab.dl.framework.domains.prediction.PredictionInterpreter;
import org.campagnelab.dl.genotype.predictions.TrueGenotypeOutputLayerPrediction;
import org.campagnelab.dl.varanalysis.protobuf.BaseInformationRecords;
import org.nd4j.linalg.api.complex.IComplexNumber;
import org.nd4j.linalg.api.ndarray.INDArray;
import org.nd4j.linalg.api.ops.impl.accum.Max;
import org.nd4j.linalg.api.ops.impl.accum.Mean;
import org.nd4j.linalg.api.ops.impl.accum.Sum;
import org.nd4j.linalg.api.ops.impl.indexaccum.IMax;
import org.nd4j.linalg.factory.Nd4j;
import org.nd4j.linalg.indexing.NDArrayIndex;

/**
 * Created by joshuacohen on 2/10/17.
 */
public class TrueGenotypeOutputLayerInterpreter implements
        PredictionInterpreter<BaseInformationRecords.BaseInformation, TrueGenotypeOutputLayerPrediction> {

    // TODO: Check more thoroughly if averaging probabilities makes sense for overallProbability
    @Override
    public TrueGenotypeOutputLayerPrediction interpret(INDArray trueLabels, INDArray output, int predictionIndex) {
        TrueGenotypeOutputLayerPrediction trueGenotypePrediction = new TrueGenotypeOutputLayerPrediction();
        INDArray trueLabelForRecord = trueLabels.getRow(predictionIndex);
        INDArray predictedLabelForRecord = output.getRow(predictionIndex);
        String trueGenotype = getGenotypeFromINDArray(trueLabelForRecord);
        String predictedGenotype = getGenotypeFromINDArray(predictedLabelForRecord);
        trueGenotypePrediction.trueGenotype = trueGenotype;
        trueGenotypePrediction.predictedGenotype = predictedGenotype;
        trueGenotypePrediction.isPredictedIndel = predictedGenotype.length() > 3;
        trueGenotypePrediction.overallProbability = getAverageFloatMaxArray(predictedLabelForRecord);;
        return trueGenotypePrediction;
    }

    @Override
    public TrueGenotypeOutputLayerPrediction interpret(BaseInformationRecords.BaseInformation record, INDArray output) {
        TrueGenotypeOutputLayerPrediction trueGenotypePrediction = new TrueGenotypeOutputLayerPrediction();
        String predictedGenotype = getGenotypeFromINDArray(output);
        trueGenotypePrediction.inspectRecord(record);
        trueGenotypePrediction.predictedGenotype = predictedGenotype;
        trueGenotypePrediction.isPredictedIndel = predictedGenotype.length() > 3;
        trueGenotypePrediction.overallProbability = getAverageFloatMaxArray(output);
        return trueGenotypePrediction;
    }

    private static String getGenotypeFromINDArray(INDArray genotypeArray) {
        INDArray genotypeLabels = getIntArgMaxArray(genotypeArray);
        StringBuilder genotypeBuilder = new StringBuilder();
        for (int i = 0; i < genotypeLabels.length(); i++) {
            int genotypeLabel = genotypeLabels.getInt(i);
            genotypeBuilder.append(labelToBase(genotypeLabel));
        }
        return genotypeBuilder.toString();
    }

    private static INDArray getIntArgMaxArray(INDArray array) {
        int maxValidIndex = Nd4j.getExecutioner().exec(new Sum(array), 0).gt(0).sumNumber().intValue();
        INDArray argMax = Nd4j.getExecutioner().exec(new IMax(array), 0);
        return maxValidIndex > 0
                ? argMax.get(NDArrayIndex.all(), NDArrayIndex.interval(0, maxValidIndex))
                : argMax;
    }

    private static double getAverageFloatMaxArray(INDArray array) {
        INDArray max = Nd4j.getExecutioner().exec(new Max(array), 0);
        int maxValidIndex = max.gt(0).sumNumber().intValue();
        INDArray truncatedMax = maxValidIndex > 0
                ? max.get(NDArrayIndex.all(), NDArrayIndex.interval(0, maxValidIndex))
                : max;
        return Nd4j.getExecutioner().execAndReturn(new Mean(truncatedMax)).getFinalResult().doubleValue();
    }

    private static char labelToBase(int label) {
        switch (label) {
            case 0:
                return 'A';
            case 1:
                return 'T';
            case 2:
                return 'C';
            case 3:
                return 'G';
            case 4:
                return 'N';
            case 5:
                return '-';
            case 6:
                return '|';
            case 7:
                return '*';
            case 8:
                return '$';
            default:
                return ' ';
        }
    }
}
