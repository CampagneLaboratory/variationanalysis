package org.campagnelab.dl.varanalysis.learning.iterators;

import it.unimi.dsi.fastutil.objects.ObjectArrayList;
import org.campagnelab.dl.varanalysis.learning.genotypes.GenotypeCountFactory;
import org.campagnelab.dl.model.utils.FeatureMapper;
import org.campagnelab.dl.varanalysis.learning.mappers.GenotypeCount;
import org.campagnelab.dl.varanalysis.protobuf.BaseInformationRecords;

import java.util.Collections;
/**
 * AbstractFeatureMapper encapsulates behavior common to many feature mappers.
 * Created by fac2003 on 6/3/16.
 *
 * @author Fabien Campagne
 */
public abstract class AbstractFeatureMapper implements FeatureMapper {
    public static final int MAX_GENOTYPES = 5;
    private static final int N_GENOTYPE_INDEX = 6;

    private boolean oneSampleHasTumor(java.util.List<org.campagnelab.dl.varanalysis.protobuf.BaseInformationRecords.SampleInfo> samples) {
        for (BaseInformationRecords.SampleInfo sample : samples) {
            if (sample.getIsTumor()) return true;
        }
        return false;

    }

    protected ObjectArrayList<? extends GenotypeCount> getAllCounts(BaseInformationRecords.BaseInformationOrBuilder record, GenotypeCountFactory factory, boolean isTumor) {
        return getAllCounts(record, factory, isTumor, true);
    }

    protected ObjectArrayList<? extends GenotypeCount> getAllCounts(BaseInformationRecords.BaseInformationOrBuilder record, GenotypeCountFactory factory, boolean isTumor, boolean sort) {
        int sampleIndex = isTumor ? 1 : 0;
        ObjectArrayList<GenotypeCount> list = new ObjectArrayList<>();
        int genotypeIndex = 0;
        for (int i = 0; i < record.getSamples(0).getCountsCount(); i++){
            int germCount = record.getSamples(0).getCounts(i).getGenotypeCountForwardStrand() + record.getSamples(0).getCounts(i).getGenotypeCountReverseStrand();
            BaseInformationRecords.CountInfo genoInfo = record.getSamples(sampleIndex).getCounts(i);
            int forwCount = genoInfo.getGenotypeCountForwardStrand();
            int revCount = genoInfo.getGenotypeCountReverseStrand();
            GenotypeCount count = factory.create();
            count.set(forwCount, revCount, genoInfo.getToSequence(),i,germCount);
            initializeCount(record.getSamples(sampleIndex).getCounts(i), count);
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

    private BaseInformationRecords.BaseInformationOrBuilder recordCached[][] = new BaseInformationRecords.BaseInformationOrBuilder[2][2];
    private ObjectArrayList<? extends GenotypeCount> cachedResult[][] = new ObjectArrayList[2][2];

    protected ObjectArrayList<? extends GenotypeCount> getAllCounts(BaseInformationRecords.BaseInformationOrBuilder record,
                                                                    boolean isTumor, boolean sort) {
        ObjectArrayList<? extends GenotypeCount> cached = getCachedResult(isTumor, sort);
        if (cached != null && record.equals(recordCached[isTumor ? 1 : 0][sort ? 1 : 0])) {
            return cached;
        } else {

            assert oneSampleHasTumor(record.getSamplesList()) : "at least one sample must have hasTumor=true.";

            for (int i = 0; i < record.getSamplesCount(); i++) {
                if (isTumor != record.getSamples(i).getIsTumor()) continue;
                // a subclass is expected to override getGenotypeCountFactory to provide its own type for Genotype counts:
                cached = getAllCounts(record, getGenotypeCountFactory(), isTumor, sort);
                putInCache(record, cached, isTumor, sort);
                return cached;
            }
            throw new InternalError("At least one sample matching isTumor, and one matching not isTumor must be found.");
        }
    }

    private void putInCache(BaseInformationRecords.BaseInformationOrBuilder record, ObjectArrayList<? extends GenotypeCount> cached, boolean isTumor, boolean sort) {
        int index1 = isTumor ? 1 : 0;
        int index2 = sort ? 1 : 0;
        recordCached[index1][index2] = record;
        cachedResult[index1][index2] = cached;
    }

    private ObjectArrayList<? extends GenotypeCount> getCachedResult(boolean isTumor, boolean sort) {
        int index1 = isTumor ? 1 : 0;
        int index2 = sort ? 1 : 0;
        return cachedResult[index1][index2];
    }

    protected ObjectArrayList<? extends GenotypeCount> getAllCounts(BaseInformationRecords.BaseInformationOrBuilder record,
                                                                    boolean isTumor) {
        return getAllCounts(record, isTumor, true);
    }

    protected abstract GenotypeCountFactory getGenotypeCountFactory();
}
