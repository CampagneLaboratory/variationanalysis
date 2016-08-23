package org.campagnelab.dl.model.utils.mappers.trio;

import it.unimi.dsi.fastutil.objects.ObjectArrayList;
import org.campagnelab.dl.model.utils.genotypes.GenotypeCountFactory;
import org.campagnelab.dl.model.utils.mappers.AbstractFeatureMapper;
import org.campagnelab.dl.model.utils.mappers.FeatureMapper;
import org.campagnelab.dl.model.utils.mappers.GenotypeCount;
import org.campagnelab.dl.varanalysis.protobuf.BaseInformationRecords;

import java.util.Collections;

/**
 * AbstractFeatureMapper encapsulates behavior common to many feature mappers.
 * Created by fac2003 on 6/3/16.
 *
 * @author Fabien Campagnei
 */
public abstract class AbstractFeatureMapperTrio extends AbstractFeatureMapper {

    protected ObjectArrayList<? extends GenotypeCount> getAllCounts(BaseInformationRecords.BaseInformationOrBuilder record, GenotypeCountFactory factory, int sampleIndex, boolean sort) {
        ObjectArrayList<GenotypeCount> list = new ObjectArrayList<>();
        int genotypeIndex = 0;
        for (int i = 0; i < record.getSamples(0).getCountsCount(); i++){
            //use somatic counts for compare
            int compareCount = record.getSamples(2).getCounts(i).getGenotypeCountForwardStrand() + record.getSamples(2).getCounts(i).getGenotypeCountReverseStrand();
            BaseInformationRecords.CountInfo genoInfo = record.getSamples(sampleIndex).getCounts(i);
            int forwCount = genoInfo.getGenotypeCountForwardStrand();
            int revCount = genoInfo.getGenotypeCountReverseStrand();
            GenotypeCount count = factory.create();
            count.set(forwCount, revCount, genoInfo.getToSequence(),i,compareCount);
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

    private BaseInformationRecords.BaseInformationOrBuilder recordCached[][] = new BaseInformationRecords.BaseInformationOrBuilder[3][2];
    private ObjectArrayList<? extends GenotypeCount> cachedResult[][] = new ObjectArrayList[3][2];

    protected ObjectArrayList<? extends GenotypeCount> getAllCounts(BaseInformationRecords.BaseInformationOrBuilder record,
                                                                    int sampleIndex, boolean sort) {
        ObjectArrayList<? extends GenotypeCount> cached = getCachedResult(sampleIndex, sort);
        if (cached != null && record.equals(recordCached[sampleIndex][sort ? 1 : 0])) {
            return cached;
        } else {

            assert oneSampleHasTumor(record.getSamplesList()) : "at least one sample must have hasTumor=true.";

            // a subclass is expected to override getGenotypeCountFactory to provide its own type for Genotype counts:
            cached = getAllCounts(record, getGenotypeCountFactory(), sampleIndex, sort);
            putInCache(record, cached, sampleIndex, sort);
            return cached;
        }
    }

    private void putInCache(BaseInformationRecords.BaseInformationOrBuilder record, ObjectArrayList<? extends GenotypeCount> cached, int sampleIndex, boolean sort) {
        int index2 = sort ? 1 : 0;
        recordCached[sampleIndex][index2] = record;
        cachedResult[sampleIndex][index2] = cached;
    }

    private ObjectArrayList<? extends GenotypeCount> getCachedResult(int sampleIndex, boolean sort) {
        int index2 = sort ? 1 : 0;
        return cachedResult[sampleIndex][index2];
    }

    protected ObjectArrayList<? extends GenotypeCount> getAllCounts(BaseInformationRecords.BaseInformationOrBuilder record,
                                                                    int sampleIndex) {
        return getAllCounts(record, sampleIndex, true);
    }

    protected abstract GenotypeCountFactory getGenotypeCountFactory();
}
