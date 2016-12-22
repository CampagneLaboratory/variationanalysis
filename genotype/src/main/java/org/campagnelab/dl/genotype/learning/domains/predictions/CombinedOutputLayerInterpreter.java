package org.campagnelab.dl.genotype.learning.domains.predictions;

import org.campagnelab.dl.framework.domains.prediction.PredictionInterpreter;
import org.campagnelab.dl.genotype.mappers.CombinedLabelsMapper;
import org.campagnelab.dl.genotype.predictions.CombinedGenotypePrediction;
import org.campagnelab.dl.genotype.predictions.CombinedOutputLayerPrediction;
import org.campagnelab.dl.varanalysis.protobuf.BaseInformationRecords;
import org.nd4j.linalg.api.ndarray.INDArray;

/**
 * Interpret model output to determine if a site is homozygous/het and what hte homosygous genotype is.
 * Created by rct66 on 11/12/16.
 */
public class CombinedOutputLayerInterpreter extends SortingCountInterpreter<CombinedOutputLayerPrediction>
        implements PredictionInterpreter<BaseInformationRecords.BaseInformation, CombinedOutputLayerPrediction> {


    public CombinedOutputLayerInterpreter() {
        super(true);
    }

    double probability;

    @Override
    public CombinedOutputLayerPrediction interpret(INDArray trueLabels, INDArray output, int predictionIndex) {
        throw new RuntimeException("a wrong interpret method was called on the homozygous interpeter");
    }

    @Override
    public CombinedOutputLayerPrediction interpret(BaseInformationRecords.BaseInformation record, INDArray output) {
        CombinedOutputLayerPrediction pred = new CombinedOutputLayerPrediction();
        pred.inspectRecord(record);
        pred.predictedGenotype = getPrediction(record, output);
        pred.overallProbability = probability;
        return pred;
    }

    public String getPrediction(BaseInformationRecords.BaseInformation currentRecord, INDArray output) {
        sortedCountRecord = sort(currentRecord);
        probability = -1;
        int maxIndex = -1;
        for (int i = 0; i < CombinedLabelsMapper.NUM_LABELS; i++) {
            double outputDouble = output.getDouble(0, i);
            if (probability < outputDouble) {
                maxIndex = i;
                probability = outputDouble;
            }
        }

        switch (maxIndex) {
            case 0:
                String homozygAllele = sortedCountRecord.getSamples(0).getCounts(0).getToSequence();
                return homozygAllele + "/" + homozygAllele;
            case 1:
                String firstAllele = sortedCountRecord.getSamples(0).getCounts(0).getToSequence();
                String secondAllele = sortedCountRecord.getSamples(0).getCounts(1).getToSequence();
                return firstAllele + "/" + secondAllele;
            case 2:
                homozygAllele = sortedCountRecord.getSamples(0).getCounts(1).getToSequence();
                return homozygAllele + "/" + homozygAllele;
            default:
                // other site predicted, mark as wrong
                return "_";
        }
    }
}
