package org.campagnelab.dl.somatic.learning.domains;

import org.campagnelab.dl.framework.domains.prediction.PredictionInterpreter;
import org.campagnelab.dl.somatic.learning.domains.predictions.IsMutatedBasePrediction;
import org.campagnelab.dl.somatic.learning.domains.predictions.IsMutatedPrediction;
import org.campagnelab.dl.somatic.learning.domains.predictions.IsSomaticMutationInterpreter;
import org.campagnelab.dl.varanalysis.protobuf.BaseInformationRecords;
import org.nd4j.linalg.api.ndarray.INDArray;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

/**
 * Interpret whether the output indicates the site is mutated. Used with IsBaseMutatedMapper.
 * Created by fac2003 on 12/5/16.
 */
public class IsBaseMutatedInterpreter extends IsSomaticMutationInterpreter implements PredictionInterpreter<BaseInformationRecords.BaseInformation, IsMutatedPrediction> {

    @Override
    public IsMutatedBasePrediction interpret(INDArray trueLabels, INDArray output, int exampleIndex) {
        IsMutatedBasePrediction prediction = new IsMutatedBasePrediction();
        prediction.trueLabelYes = 1 - trueLabels.getDouble(exampleIndex, 0);

        double probabilityNotMutated = output.getDouble(exampleIndex, 0);
        prediction.predictedLabelNo = probabilityNotMutated;
        prediction.predictedLabelYes = 1 - probabilityNotMutated;

        // represent the allele as 0, 1, according to index in sorted counts:
        prediction.predictedMutatedAllele = Integer.toString(getArgMaxIndex(output, exampleIndex) - 1);
        prediction.trueMutatedAllele = Integer.toString(getArgMaxIndex(trueLabels, exampleIndex) - 1);

        prediction.predictedLabelNo = probabilityNotMutated;
        prediction.predictedLabelYes = 1 - probabilityNotMutated;

        return prediction;
    }

    private int getArgMaxIndex(INDArray output, int predictionIndex) {
        double max = Double.NEGATIVE_INFINITY;
        int argMaxIndex = -1;
        for (int i = 0; i < output.size(1); i++) {
            double value = output.getDouble(predictionIndex, i);
            if (value > max) {
                max = value;
                argMaxIndex = i;
            }
        }
        return argMaxIndex;
    }

    @Override
    public IsMutatedBasePrediction interpret(BaseInformationRecords.BaseInformation record, INDArray output) {
        IsMutatedBasePrediction prediction = new IsMutatedBasePrediction();
        prediction.inspectRecord(record);
        int predictionIndex = 0;

        prediction.trueLabelYes = record.getMutated() ? 1.0 : 0.0;
        double probabilityNotMutated = output.getDouble(predictionIndex, 0);
        prediction.predictedLabelNo = probabilityNotMutated;
        prediction.predictedLabelYes = 1 - probabilityNotMutated;

        // represent the allele as 0, 1, according to index in sorted counts:
        prediction.predictedLabelNo = probabilityNotMutated;
        prediction.predictedLabelYes = 1 - probabilityNotMutated;

        int argMaxIndex = getArgMaxIndex(output, predictionIndex);

        assert argMaxIndex != -1 : "argMaxIndex must be set";
        prediction.sortedBaseIndex = argMaxIndex - 1;
        // determine the base allele for called somatic allele:
        final List<BaseInformationRecords.CountInfo> counts = record.getSamples(record.getSamplesCount() - 1).getCountsList();
        ArrayList<BaseInformationRecords.CountInfo> sorted = new ArrayList<>();
        sorted.addAll(counts);
        Collections.sort(sorted, (o1, o2) ->
                (o2.getGenotypeCountForwardStrand() + o2.getGenotypeCountReverseStrand()) - (o1.getGenotypeCountForwardStrand() + o1.getGenotypeCountReverseStrand())
        );
        if (argMaxIndex != 0) {
            prediction.predictedMutatedAllele = sorted.get(argMaxIndex - 1).getToSequence();
        } else {
            prediction.predictedMutatedAllele = ".";
        }
        prediction.trueMutatedAllele = record.getMutatedBase();
        return prediction;
    }
}
