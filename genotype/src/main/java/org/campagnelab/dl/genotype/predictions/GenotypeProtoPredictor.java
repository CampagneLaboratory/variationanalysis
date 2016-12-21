package org.campagnelab.dl.genotype.predictions;

import org.campagnelab.dl.framework.domains.DomainDescriptor;
import org.campagnelab.dl.framework.mappers.FeatureMapper;
import org.campagnelab.dl.framework.models.ModelOutputHelper;
import org.campagnelab.dl.genotype.learning.domains.predictions.HomozygousInterpreter;
import org.campagnelab.dl.genotype.learning.domains.predictions.HomozygousPrediction;
import org.campagnelab.dl.genotype.learning.domains.predictions.SingleGenotypeInterpreter;
import org.campagnelab.dl.genotype.learning.domains.predictions.SingleGenotypePrediction;
import org.campagnelab.dl.varanalysis.protobuf.BaseInformationRecords;
import org.deeplearning4j.nn.api.Model;
import org.nd4j.linalg.api.ndarray.INDArray;

/**
 * Created by fac2003 on 12/18/16.
 */
public class GenotypeProtoPredictor {

    private static int HOMOZYGOUS_OUTPUT_INDEX;
    private final DomainDescriptor domainDescriptor;
    private final Model model;
    private final FeatureMapper mapper;
    private HomozygousInterpreter isHomozygous;
    private SingleGenotypeInterpreter[] singleGenotypeInterpreters;
    private int[] singleGenotypeOutputIndices;
    private ModelOutputHelper outputHelper;

    public GenotypeProtoPredictor(DomainDescriptor domainDescriptor, Model model, FeatureMapper featureMapper) {
        this.domainDescriptor = domainDescriptor;
        this.model = model;
        this.mapper = featureMapper;
        this.outputHelper = new ModelOutputHelper();
        String[] outputNames = domainDescriptor.getComputationalGraph().getOutputNames();

        singleGenotypeOutputIndices = new int[outputNames.length - 1];
        singleGenotypeInterpreters = new SingleGenotypeInterpreter[outputNames.length - 1];
        int outputIndex = 0;
        int singleGenotypeOutputIndex = 0;
        for (String outputName : outputNames) {
            switch (outputName) {
                case "homozygous":
                    isHomozygous = (HomozygousInterpreter) domainDescriptor.getPredictionInterpreter("homozygous");
                    HOMOZYGOUS_OUTPUT_INDEX = outputIndex++;
                    break;
                default:
                    singleGenotypeOutputIndices[singleGenotypeOutputIndex] = outputIndex;
                    singleGenotypeInterpreters[singleGenotypeOutputIndex++] = (SingleGenotypeInterpreter) domainDescriptor.getPredictionInterpreter(outputName);
                    break;
            }
        }
    }

    public GenotypePrediction5Out predictGenotype(BaseInformationRecords.BaseInformation record) {
        assert model != null : "Model cannot be null";
        INDArray arrayPredicted = null;
        GenotypePrediction5Out prediction = new GenotypePrediction5Out();

        outputHelper.predictForNextRecord(model, record, mapper);

        HomozygousPrediction homozygousPrediction = isHomozygous.interpret(record, outputHelper.getOutput(HOMOZYGOUS_OUTPUT_INDEX));
        SingleGenotypePrediction[] singleGenotypePredictions = new SingleGenotypePrediction[singleGenotypeOutputIndices.length];
        int i = 0;
        for (int singleGenotypeIndex : singleGenotypeOutputIndices) {
            singleGenotypePredictions[i] = singleGenotypeInterpreters[i].interpret(record, outputHelper.getOutput(singleGenotypeIndex));
            i++;
        }
        prediction.set(homozygousPrediction, singleGenotypePredictions);
        return prediction;
    }
}
