package org.campagnelab.dl.genotype.predictions;

import org.campagnelab.dl.framework.domains.DomainDescriptor;
import org.campagnelab.dl.framework.domains.prediction.Prediction;
import org.campagnelab.dl.framework.domains.prediction.PredictionInterpreter;
import org.campagnelab.dl.framework.mappers.FeatureMapper;
import org.campagnelab.dl.framework.models.ModelOutputHelper;
import org.campagnelab.dl.genotype.learning.domains.predictions.HomozygousInterpreter;
import org.campagnelab.dl.genotype.learning.domains.predictions.SingleGenotypeInterpreter;
import org.campagnelab.dl.varanalysis.protobuf.BaseInformationRecords;
import org.deeplearning4j.nn.api.Model;
import org.nd4j.linalg.api.ndarray.INDArray;

import java.util.ArrayList;
import java.util.List;

/**
 * Created by fac2003 on 12/18/16.
 */
public class GenotypeProtoPredictor {

    private final DomainDescriptor domainDescriptor;
    private final Model model;
    private final FeatureMapper mapper;
    private ModelOutputHelper outputHelper;

    public GenotypeProtoPredictor(DomainDescriptor domainDescriptor, Model model, FeatureMapper featureMapper) {
        this.domainDescriptor = domainDescriptor;
        this.model = model;
        this.mapper = featureMapper;
        this.outputHelper = new ModelOutputHelper();
        String[] outputNames = domainDescriptor.getComputationalGraph().getOutputNames();
        int outputIndex = 0;

        interpretors = new PredictionInterpreter[outputNames.length];
        for (String outputName : outputNames) {
            interpretors[outputIndex++] = domainDescriptor.getPredictionInterpreter(outputName);
        }


    }

    private PredictionInterpreter[] interpretors;
    private List<Prediction> predictions = new ArrayList<>();

    public GenotypePrediction predictGenotype(BaseInformationRecords.BaseInformation currentRecord) {
        assert model != null : "Model cannot be null";

        outputHelper.predictForNextRecord(model, currentRecord, mapper);
        predictions.clear();
        for (int outputIndex = 0; outputIndex < domainDescriptor.getNumModelOutputs(); outputIndex++) {
            INDArray outputPredictions = outputHelper.getOutput(outputIndex);

            if (interpretors[outputIndex] != null) {
                Prediction prediction = interpretors[outputIndex].interpret(currentRecord, outputPredictions);
                prediction.outputIndex = outputIndex;
                predictions.add(prediction);
            }
        }

        GenotypePrediction overallPrediction = (GenotypePrediction) domainDescriptor.aggregatePredictions(predictions);

        return overallPrediction;
    }
}
