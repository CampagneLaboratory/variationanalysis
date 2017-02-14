package org.campagnelab.dl.genotype.predictions;

import org.campagnelab.dl.framework.domains.prediction.PredictionInterpreter;
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
public class TrueGenotypeInterpreter implements
        PredictionInterpreter<BaseInformationRecords.BaseInformation, TrueGenotypePrediction> {

    // TODO: In methods below, check how predictedGenotype should be set
    // TODO: In methods below, check how isVariantProbability should be set- both currently set to overallProbability
    // TODO: Check more thoroughly if averaging probabilities makes sense

    // TODO: In below method, seemingly no way to check trueFrom or predictedFrom- just set to genotype
    @Override
    public TrueGenotypePrediction interpret(INDArray trueLabels, INDArray output, int predictionIndex) {
        TrueGenotypePrediction trueGenotypePrediction = new TrueGenotypePrediction();
        INDArray trueLabelForRecord = trueLabels.getRow(predictionIndex);
        INDArray predictedLabelForRecord = output.getRow(predictionIndex);
        String trueGenotype = getGenotypeFromINDArray(trueLabelForRecord);
        String predictedGenotype = getGenotypeFromINDArray(predictedLabelForRecord);
        trueGenotypePrediction.trueGenotype = trueGenotype;
        trueGenotypePrediction.isIndel = trueGenotype.length() > 1;
        trueGenotypePrediction.predictedGenotype = predictedGenotype;
        trueGenotypePrediction.isPredictedIndel = predictedGenotype.length() > 1;
        trueGenotypePrediction.trueFrom = trueGenotype;
        trueGenotypePrediction.predictedFrom = predictedGenotype;
        double predictionProbability = getAverageFloatMaxArray(predictedLabelForRecord);
        trueGenotypePrediction.overallProbability = predictionProbability;
        trueGenotypePrediction.isVariantProbability = predictionProbability;
        return trueGenotypePrediction;
    }

    // TODO: In below method, predictedFrom is just set using the default from inspectRecord
    @Override
    public TrueGenotypePrediction interpret(BaseInformationRecords.BaseInformation record, INDArray output) {
        TrueGenotypePrediction trueGenotypePrediction = new TrueGenotypePrediction();
        trueGenotypePrediction.inspectRecord(record);
        String predictedGenotype = getGenotypeFromINDArray(output);
        trueGenotypePrediction.predictedGenotype = predictedGenotype;
        trueGenotypePrediction.isPredictedIndel = predictedGenotype.length() > 1;
        double predictionProbability = getAverageFloatMaxArray(output);
        trueGenotypePrediction.overallProbability = predictionProbability;
        trueGenotypePrediction.isVariantProbability = predictionProbability;
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
        return argMax.get(NDArrayIndex.all(), NDArrayIndex.interval(0, maxValidIndex));
    }

    private static  double getAverageFloatMaxArray(INDArray array) {
        INDArray max = Nd4j.getExecutioner().exec(new Max(array), 0);
        return Nd4j.getExecutioner().execAndReturn(new Mean(max)).getFinalResult().doubleValue();
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
                return '/';
            default:
                return ' ';
        }
    }
}
