package org.campagnelab.dl.genotype.segments;

import it.unimi.dsi.fastutil.floats.FloatArrayList;
import it.unimi.dsi.fastutil.floats.FloatList;
import org.campagnelab.dl.framework.mappers.FeatureMapper;
import org.campagnelab.dl.genotype.helpers.GenotypeHelper;
import org.campagnelab.dl.genotype.mappers.SingleBaseMapper;
import org.campagnelab.dl.genotype.tools.SBIToSSIConverterArguments;
import org.campagnelab.dl.varanalysis.protobuf.BaseInformationRecords;
import org.campagnelab.dl.varanalysis.protobuf.SegmentInformationRecords;

/**
 * Function used by the sbi to ssi converter to populate features and labels.
 */
public class MyFillInFeaturesFunction implements FillInFeaturesFunction {
    private final FeatureMapper featureMapper;
    private final SegmentLabelMapper labelMapper;
    private final SBIToSSIConverterArguments arguments;

    public MyFillInFeaturesFunction(FeatureMapper featureMapper, SegmentLabelMapper labelMapper, SBIToSSIConverterArguments arguments) {
        this.featureMapper = featureMapper;
        this.labelMapper = labelMapper;
        this.arguments = arguments;
    }

    // any genomic site that has strictly more indel supporting reads than the below threshold will be marked has candidateIndel.
    private int candidateIndelThreshold = 0;

    @Override
    public SegmentInformationRecords.Base.Builder apply(BaseInformationRecords.BaseInformation baseInformation) {

        SegmentInformationRecords.Base.Builder builder = SegmentInformationRecords.Base.newBuilder();
        String trueGenotype = baseInformation.getTrueGenotype();
        builder.addAllTrueLabel(GenotypeHelper.getAlleles(trueGenotype));

        builder.setHasCandidateIndel(hasCandidateIndel(baseInformation));
        builder.setHasTrueIndel(
                GenotypeHelper.isIndel(baseInformation.getReferenceBase(), baseInformation.getTrueGenotype()));
        builder.setIsVariant(
                GenotypeHelper.isVariant(baseInformation.getTrueGenotype(), baseInformation.getReferenceBase()));
        builder.setReferenceAllele(baseInformation.getReferenceBase());
        builder.setFormattedCounts(FormatterCountHelper.format(baseInformation.getSamples(0)));
        if (args().mapFeatures) {
            FloatList features = new FloatArrayList(featureMapper.numberOfFeatures());

            featureMapper.prepareToNormalize(baseInformation, 0);
            if (trueGenotype.length() > 3) {
                //    System.out.println("Indel:" + baseInformation.getTrueGenotype());
            }
            features.clear();
            for (int featureIndex = 0; featureIndex < featureMapper.numberOfFeatures(); featureIndex++) {
                features.add(featureMapper.produceFeature(baseInformation, featureIndex));
            }
            builder.clearFeatures();
            builder.addAllFeatures(features);
        }
        if (args().mapLabels) {
            if (trueGenotype.length() == 1) {
                trueGenotype = trueGenotype + "|" + trueGenotype;
            }
            if (trueGenotype.length() > 3) {
                trueGenotype = trimTrueGenotype(trueGenotype);
            }
            float[] labels = labelMapper.map(trueGenotype.replaceAll("\\|", "/"));
            builder.clearLabels();
            for (float labelValue : labels) {
                builder.addLabels(labelValue);
            }
        }
        return builder;
    }

    private SBIToSSIConverterArguments args() {
        return arguments;
    }

    private boolean hasCandidateIndel(BaseInformationRecords.BaseInformation baseInformation) {
        return SegmentUtil.hasCandidateIndel(baseInformation, candidateIndelThreshold);
    }

    private String trimTrueGenotype(String trueGenotype) {
        int secondBaseIndex = trueGenotype.indexOf("|");
        final String trimmed = trueGenotype.charAt(0) + "|" + trueGenotype.charAt(secondBaseIndex - 1);
        return trimmed;
    }
}
