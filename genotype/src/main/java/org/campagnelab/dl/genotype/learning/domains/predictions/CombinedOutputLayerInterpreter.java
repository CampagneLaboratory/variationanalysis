package org.campagnelab.dl.genotype.learning.domains.predictions;

import org.campagnelab.dl.framework.domains.prediction.PredictionInterpreter;
import org.campagnelab.dl.genotype.mappers.CombinedLabelsMapper;
import org.campagnelab.dl.genotype.predictions.CombinedOutputLayerPrediction;
import org.campagnelab.dl.varanalysis.protobuf.BaseInformationRecords;
import org.jetbrains.annotations.NotNull;
import org.nd4j.linalg.api.ndarray.INDArray;

/**
 * Interpret model output to determine if a site is homozygous/het and what hte homosygous genotype is.
 * Created by rct66 on 11/12/16.
 */
public class CombinedOutputLayerInterpreter extends SortingCountInterpreter<CombinedOutputLayerPrediction>
        implements PredictionInterpreter<BaseInformationRecords.BaseInformation, CombinedOutputLayerPrediction> {

    private static final int MAX_GENOTYPES = CombinedLabelsMapper.NUM_LABELS;
    private String[] toSequences;
    private double probability;

    public CombinedOutputLayerInterpreter() {
        super(true);
        // the default toSequences are integers, because we don't have access to exact genotypes (we have no record),
        // but integers are sufficient to calculate performance stats.
        toSequences = new String[MAX_GENOTYPES];
        for (int genotypeIndex = 0; genotypeIndex < MAX_GENOTYPES; genotypeIndex++) {
            toSequences[genotypeIndex] = Integer.toString(genotypeIndex);
        }
        // the other genotype can never match:
        toSequences[3]="-1";
    }


    @Override
    public CombinedOutputLayerPrediction interpret(INDArray trueLabels, INDArray output, int predictionIndex) {
        CombinedOutputLayerPrediction p = new CombinedOutputLayerPrediction();
        p.trueGenotype = reconstructGenotype(trueLabels);
        p.predictedGenotype = reconstructGenotype(output);
        return p;
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
        for (int genotypeIndex = 0; genotypeIndex < MAX_GENOTYPES; genotypeIndex++) {
            toSequences[genotypeIndex] = sortedCountRecord.getSamples(0).getCounts(genotypeIndex).getToSequence();
        }
        return reconstructGenotype(output);
    }

    @NotNull
    private String reconstructGenotype(INDArray output) {
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
                String homozygAllele = toSequences[0];
                return homozygAllele + "/" + homozygAllele;
            case 1:
                String firstAllele = toSequences[0];
                String secondAllele = toSequences[1];
                return firstAllele + "/" + secondAllele;
            case 2:
                homozygAllele = toSequences[1];
                return homozygAllele + "/" + homozygAllele;
            default:
                // other site predicted, mark as wrong
                return "_";
        }
    }
}
