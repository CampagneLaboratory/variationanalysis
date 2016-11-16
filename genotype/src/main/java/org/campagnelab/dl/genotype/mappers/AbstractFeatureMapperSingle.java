package org.campagnelab.dl.genotype.mappers;

import it.unimi.dsi.fastutil.objects.ObjectArrayList;
import org.campagnelab.dl.framework.mappers.FeatureNameMapper;
import org.campagnelab.dl.somatic.genotypes.GenotypeCountFactory;
import org.campagnelab.dl.somatic.mappers.GenotypeCount;
import org.campagnelab.dl.varanalysis.protobuf.BaseInformationRecords;
import org.nd4j.linalg.api.ndarray.INDArray;

import java.util.Arrays;
import java.util.Collections;

/**
 * AbstractFeatureMapper encapsulates behavior common to many feature mappers.
 * Created by fac2003 on 6/3/16.
 *
 * @author Fabien Campagne
 */
public abstract class AbstractFeatureMapperSingle<T extends BaseInformationRecords.BaseInformationOrBuilder > implements FeatureNameMapper<T> {
    public static final int MAX_GENOTYPES = 5;
    public static final int N_GENOTYPE_INDEX = 6;
    private float[] buffer;
    private int sampleIndex;

    protected float[] getBuffer() {

        if (buffer==null) {
            buffer = new float[numberOfFeatures()];
        }else{
            Arrays.fill(buffer,0f);
        }
        return buffer;
    }
    @Override
    public String getFeatureName(int featureIndex) {
        return null;
    }

    protected ObjectArrayList<? extends GenotypeCount> getAllCounts(BaseInformationRecords.BaseInformationOrBuilder record, GenotypeCountFactory factory) {
        return getAllCounts(record, factory, true);
    }

    protected ObjectArrayList<? extends GenotypeCount> getAllCounts(BaseInformationRecords.BaseInformationOrBuilder record, GenotypeCountFactory factory, boolean sort) {
        ObjectArrayList<GenotypeCount> list = new ObjectArrayList<>();
        int genotypeIndex = 0;
        BaseInformationRecords.SampleInfo sample = record.getSamples(sampleIndex);
        for (int i = 0; i < record.getSamples(sampleIndex).getCountsCount(); i++){
            BaseInformationRecords.CountInfo genotype = sample.getCounts(i);
            int forwCount = genotype.getGenotypeCountForwardStrand();
            int revCount = genotype.getGenotypeCountReverseStrand();
            int sumCount = forwCount + revCount;
            GenotypeCount count = factory.create();
            count.set(forwCount, revCount, genotype.getToSequence(),i,sumCount);
            initializeCount(genotype, count);
            list.add(count);
        }
        // DO not increment genotypeIndex. It must remain constant for all N bases
        int genotypeIndexFor_Ns = N_GENOTYPE_INDEX;
        // pad with zero until we have 10 elements:
        while (list.size() < MAX_GENOTYPES) {
            final GenotypeCount genotypeCount = getGenotypeCountFactory().create();
            genotypeCount.set(0, 0, "N", genotypeIndexFor_Ns,0);
            list.add(genotypeCount);

        }
        //sort in decreasing order of counts:
        if (sort) {
            Collections.sort(list);
        }
        // trim the list at 5 elements because we will consider only the 5 genotypes with largest total counts:
        list.trim(MAX_GENOTYPES);

        return list;
    }


    protected abstract void initializeCount(BaseInformationRecords.CountInfo sampleCounts, GenotypeCount count);

    private BaseInformationRecords.BaseInformationOrBuilder recordCached[] = new BaseInformationRecords.BaseInformationOrBuilder[2];
    private ObjectArrayList<? extends GenotypeCount> cachedResult[] = new ObjectArrayList[2];

    protected ObjectArrayList<? extends GenotypeCount> getAllCounts(BaseInformationRecords.BaseInformationOrBuilder record, boolean sort) {
        ObjectArrayList<? extends GenotypeCount> cached = getCachedResult(sort);
        if (cached != null && record.equals(recordCached[sort ? 1 : 0])) {
            return cached;
        } else {
            cached = getAllCounts(record, getGenotypeCountFactory(), sort);
            putInCache(record, cached, sort);
            return cached;
        }
    }

    private void putInCache(BaseInformationRecords.BaseInformationOrBuilder record, ObjectArrayList<? extends GenotypeCount> cached, boolean sort) {
        int index = sort ? 1 : 0;
        recordCached[index] = record;
        cachedResult[index] = cached;
    }

    private ObjectArrayList<? extends GenotypeCount> getCachedResult(boolean sort) {
        int index = sort ? 1 : 0;
        return cachedResult[index];
    }

    protected ObjectArrayList<? extends GenotypeCount> getAllCounts(BaseInformationRecords.BaseInformationOrBuilder record) {
        return getAllCounts(record, true);
    }

    protected abstract GenotypeCountFactory getGenotypeCountFactory();

    @Override
    public boolean hasMask() {
        return false;
    }

    @Override
    public void maskFeatures(T record, INDArray mask, int indexOfRecord) {
        // do nothing implementation.
    }

    @Override
    public boolean isMasked(T record, int featureIndex) {
        return false;
    }
}
