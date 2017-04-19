package org.campagnelab.dl.genotype.predictions;

import org.campagnelab.dl.framework.domains.prediction.PredictionInterpreter;
import org.campagnelab.dl.genotype.helpers.GenotypeHelper;
import org.campagnelab.dl.genotype.mappers.MetaDataLabelMapper;
import org.campagnelab.dl.genotype.mappers.RecordCountSortHelper;
import org.campagnelab.dl.varanalysis.protobuf.BaseInformationRecords;
import org.nd4j.linalg.api.ndarray.INDArray;

/**
 * This interpreter extracts meta-data (isVariant, isIndel, etc.) from the meta-data label.
 * Created by fac2003 on 12/22/16.
 */
public class MetaDataInterpreter implements PredictionInterpreter<BaseInformationRecords.BaseInformation,
        MetadataPrediction> {

    private RecordCountSortHelper sortHelper = new RecordCountSortHelper();

    @Override
    public MetadataPrediction interpret(INDArray trueLabels, INDArray output, int predictionIndex) {
        MetadataPrediction p = new MetadataPrediction();
        p.isIndel = trueLabels.getDouble(predictionIndex, MetaDataLabelMapper.IS_INDEL_FEATURE_INDEX) == 1;
        p.isVariant = trueLabels.getDouble(predictionIndex, MetaDataLabelMapper.IS_VARIANT_FEATURE_INDEX) == 1;
        p.referenceGobyIndex = (int) trueLabels.getDouble(predictionIndex, MetaDataLabelMapper.IS_MATCHING_REF_FEATURE_INDEX);
        p.sorted2OriginalCountIndices = new int[]{
                (int) trueLabels.getDouble(predictionIndex, MetaDataLabelMapper.IS_COUNT1_ORIGINAL_INDEX_FEATURE_INDEX),
                (int) trueLabels.getDouble(predictionIndex, MetaDataLabelMapper.IS_COUNT2_ORIGINAL_INDEX_FEATURE_INDEX),
                (int) trueLabels.getDouble(predictionIndex, MetaDataLabelMapper.IS_COUNT3_ORIGINAL_INDEX_FEATURE_INDEX),
                (int) trueLabels.getDouble(predictionIndex, MetaDataLabelMapper.IS_COUNT4_ORIGINAL_INDEX_FEATURE_INDEX),
                (int) trueLabels.getDouble(predictionIndex, MetaDataLabelMapper.IS_COUNT5_ORIGINAL_INDEX_FEATURE_INDEX),
                (int) trueLabels.getDouble(predictionIndex, MetaDataLabelMapper.IS_COUNT6_ORIGINAL_INDEX_FEATURE_INDEX),
                (int) trueLabels.getDouble(predictionIndex, MetaDataLabelMapper.IS_COUNT7_ORIGINAL_INDEX_FEATURE_INDEX),
                (int) trueLabels.getDouble(predictionIndex, MetaDataLabelMapper.IS_COUNT8_ORIGINAL_INDEX_FEATURE_INDEX),
                (int) trueLabels.getDouble(predictionIndex, MetaDataLabelMapper.IS_COUNT9_ORIGINAL_INDEX_FEATURE_INDEX),
                (int) trueLabels.getDouble(predictionIndex, MetaDataLabelMapper.IS_COUNT10_ORIGINAL_INDEX_FEATURE_INDEX)

        };
        return p;
    }

    @Override
    public MetadataPrediction interpret(BaseInformationRecords.BaseInformation record, INDArray output) {
        MetadataPrediction p = new MetadataPrediction();
        BaseInformationRecords.SampleInfo sample = record.getSamples(0);
        p.isVariant = sample.getIsVariant();
        final String trueGenotype = record.getTrueGenotype();
        p.isIndel = GenotypeHelper.isIndel(record.getReferenceBase(), trueGenotype);
        p.referenceGobyIndex = MetaDataLabelMapper.calculateReferenceIndex(record);
        BaseInformationRecords.BaseInformation sortedRecord = sortHelper.sort(0, record);
        // obtain original indices for sorted counts:
        BaseInformationRecords.SampleInfo sortedCountSample = sortedRecord.getSamples(0);
        p.sorted2OriginalCountIndices = new int[]{
                sortedCountSample.getCountsCount() >= 1 ? sortedCountSample.getCounts(0).getGobyGenotypeIndex() : -1,
                sortedCountSample.getCountsCount() >= 2 ? sortedCountSample.getCounts(1).getGobyGenotypeIndex() : -1,
                sortedCountSample.getCountsCount() >= 3 ? sortedCountSample.getCounts(2).getGobyGenotypeIndex() : -1,
                sortedCountSample.getCountsCount() >= 4 ? sortedCountSample.getCounts(3).getGobyGenotypeIndex() : -1,
                sortedCountSample.getCountsCount() >= 5 ? sortedCountSample.getCounts(4).getGobyGenotypeIndex() : -1,
                sortedCountSample.getCountsCount() >= 6 ? sortedCountSample.getCounts(5).getGobyGenotypeIndex() : -1,
                sortedCountSample.getCountsCount() >= 7 ? sortedCountSample.getCounts(6).getGobyGenotypeIndex() : -1,
                sortedCountSample.getCountsCount() >= 8 ? sortedCountSample.getCounts(7).getGobyGenotypeIndex() : -1,
                sortedCountSample.getCountsCount() >= 9 ? sortedCountSample.getCounts(8).getGobyGenotypeIndex() : -1,
                sortedCountSample.getCountsCount() >= 10 ? sortedCountSample.getCounts(9).getGobyGenotypeIndex() : -1,
                sortedCountSample.getCountsCount() >= 11 ? sortedCountSample.getCounts(10).getGobyGenotypeIndex() : -1,
        };
        return p;
    }
}
