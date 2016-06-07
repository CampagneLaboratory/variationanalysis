package org.campagnelab.dl.varanalysis.learning.iterators;

import it.unimi.dsi.fastutil.objects.ObjectArrayList;
import org.campagnelab.dl.varanalysis.learning.genotypes.GenotypeCountFactory;
import org.campagnelab.dl.varanalysis.learning.mappers.FeatureMapper;
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

    private boolean oneSampleHasTumor(java.util.List<org.campagnelab.dl.varanalysis.protobuf.BaseInformationRecords.SampleInfo> samples) {
        for (BaseInformationRecords.SampleInfo sample : samples) {
            if (sample.getIsTumor()) return true;
        }
        return false;

    }

    protected ObjectArrayList<? extends GenotypeCount> getAllCounts(BaseInformationRecords.SampleInfo sample, GenotypeCountFactory factory) {
        ObjectArrayList<GenotypeCount> list = new ObjectArrayList<>();

        for (BaseInformationRecords.CountInfo sampleCounts : sample.getCountsList()) {
            GenotypeCount count = factory.create();
            count.set(sampleCounts.getGenotypeCountForwardStrand(), sampleCounts.getGenotypeCountReverseStrand(), sampleCounts.getToSequence());
            initializeCount(sampleCounts, count);
            list.add(count);
        }

        // pad with zero until we have 10 elements:
        while (list.size() < MAX_GENOTYPES) {
            final GenotypeCount genotypeCount = getGenotypeCountFactory().create();
            genotypeCount.set(0, 0, "N");
            list.add(genotypeCount);
        }
        // trim the list at 5 elements because we will consider only the 5 genotypes with largest total counts:
        list.trim(MAX_GENOTYPES);
        //sort in decreasing order of counts:
        Collections.sort(list);
        return list;
    }

    protected abstract void initializeCount(BaseInformationRecords.CountInfo sampleCounts, GenotypeCount count);

    protected ObjectArrayList<? extends GenotypeCount> getAllCounts(BaseInformationRecords.BaseInformationOrBuilder record, boolean isTumor) {
        assert oneSampleHasTumor(record.getSamplesList()) : "at least one sample must have hasTumor=true.";
        for (BaseInformationRecords.SampleInfo sampleInfo : record.getSamplesList()) {
            if (isTumor != sampleInfo.getIsTumor()) continue;
            // a subclass is expected to override getGenotypeCountFactory to provide its own type for Genotype counts:
            return getAllCounts(sampleInfo, getGenotypeCountFactory());
        }
        throw new InternalError("At least one sample matching isTumor, and one matching not isTumor must be found.");
    }

    protected abstract GenotypeCountFactory getGenotypeCountFactory();
}
