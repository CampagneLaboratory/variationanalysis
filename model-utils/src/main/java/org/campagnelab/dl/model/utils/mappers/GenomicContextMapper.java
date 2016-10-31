package org.campagnelab.dl.model.utils.mappers;

import it.unimi.dsi.fastutil.objects.ObjectArrayList;
import org.campagnelab.dl.model.utils.ProtoPredictor;
import org.campagnelab.dl.model.utils.genotypes.BaseGenotypeCountFactory;
import org.campagnelab.dl.model.utils.genotypes.GenotypeCountFactory;
import org.campagnelab.dl.model.utils.mappers.functional.TraversalHelper;
import org.campagnelab.dl.varanalysis.protobuf.BaseInformationRecords;
import org.nd4j.linalg.api.ndarray.INDArray;

import java.util.List;
import java.util.Properties;

/**
 * Maps the full genomic context using multiple onehotfeaturemapper
 * Created by rct66 on 10/25/16.
 */



public class GenomicContextMapper extends AbstractFeatureMapper implements FeatureMapper, EfficientFeatureMapper {
    private ConcatFeatureMapper delegate;
    private int contextSize;
    public GenomicContextMapper(Properties sbiProperties) {
        if (!sbiProperties.containsKey("contextSize")) {
            throw new UnsupportedOperationException("The sbip file does not contain context size information (attempted to access field 'contextSize')");
        }
        this.contextSize = Integer.parseInt(sbiProperties.getProperty("contextSize"));
        OneHotBaseMapper[] refContext = new OneHotBaseMapper[contextSize];

        for (int i = 0; i < contextSize; i++){
            refContext[i] = new OneHotBaseMapper(i, BaseInformationRecords.BaseInformationOrBuilder::getGenomicSequenceContext);
        }
        delegate = new ConcatFeatureMapper(refContext);
    }

    @Override
    public int numberOfFeatures() {
        return delegate.numberOfFeatures();

    }

    @Override
    public void prepareToNormalize(BaseInformationRecords.BaseInformationOrBuilder record, int indexOfRecord) {
        delegate.prepareToNormalize(record,indexOfRecord);
    }

    int[] indices = new int[]{0, 0};

    @Override
    public void mapFeatures(BaseInformationRecords.BaseInformationOrBuilder record, INDArray inputs, int indexOfRecord) {
        delegate.mapFeatures(record, inputs, indexOfRecord);
    }

    @Override
    public void mapFeatures(BaseInformationRecords.BaseInformationOrBuilder record, float[] inputs, int offset, int indexOfRecord) {
        for (int featureIndex = 0; featureIndex < numberOfFeatures(); featureIndex++) {
            inputs[featureIndex+offset] = produceFeature(record, featureIndex);
        }
    }

    @Override
    public float produceFeature(BaseInformationRecords.BaseInformationOrBuilder record, int featureIndex) {
        return delegate.produceFeature(record, featureIndex);
    }

    @Override
    protected void initializeCount(BaseInformationRecords.CountInfo sampleCounts, GenotypeCount count) {
        ReadIndexWithCounts myCounts = (ReadIndexWithCounts) count;
        myCounts.set(ProtoPredictor.expandFreq(sampleCounts.getReadIndicesForwardStrandList()),
                ProtoPredictor.expandFreq(sampleCounts.getReadIndicesReverseStrandList()));
    }

    @Override
    protected GenotypeCountFactory getGenotypeCountFactory() {

        return new BaseGenotypeCountFactory() {
            @Override
            public GenotypeCount create() {
                return new ReadIndexWithCounts();
            }
        };
    }

}
