package org.campagnelab.dl.genotype.mappers;

import org.campagnelab.dl.framework.mappers.ConfigurableFeatureMapper;
import org.campagnelab.dl.varanalysis.protobuf.BaseInformationRecords;
import org.campagnelab.goby.util.WarningCounter;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.util.List;
import java.util.Properties;

/**
 * Label that encodes  the combination of called alleles in a sample using one-hot encoding and a softmax for decoding.<p/>
 * <p>
 * Created by fac2003 on 2/7/17.
 */
public class SoftmaxLabelMapper extends CountSortingLabelMapper implements ConfigurableFeatureMapper {

    private WarningCounter unableToRepresent = new WarningCounter(10);
    private final float epsilon;

    public void setMaxCalledAlleles(int maxCalledAlleles) {
        this.maxCalledAlleles = maxCalledAlleles;
        MAX_VALUE = (int) Math.pow(2, maxCalledAlleles);
    }

    public int maxCalledAlleles;
    static private Logger LOG = LoggerFactory.getLogger(SoftmaxLabelMapper.class);
    private int MAX_VALUE;

    /**
     * @param sortCounts
     * @param maxCalledAlleles
     * @param epsilon          amount of label smoothing to apply.
     */
    public SoftmaxLabelMapper(int sampleIndex, boolean sortCounts, int maxCalledAlleles, float epsilon) {

        super(sortCounts);
        if (!sortCounts) LOG.warn("You should only useSoftmaxLabelMapper with unsorted counts in tests. ");

        this.epsilon = epsilon;
        setMaxCalledAlleles(maxCalledAlleles);
    }

    public SoftmaxLabelMapper(int sampleIndex, boolean sortCounts, int ploidy) {

        this(sampleIndex, sortCounts, ploidy, 0);
    }

    @Override
    public int numberOfLabels() {
        return (MAX_VALUE + 1);
    }

    private int cachedValue;

    protected void setCachedValue(boolean... isCalled) {
        cachedValue = 0;
        int index = 0;
        for (boolean called : isCalled) {
            cachedValue |= (called ? 1 : 0) << index;
            index++;
            if (cachedValue > MAX_VALUE) {
                break;
            }
        }
        if (cachedValue > MAX_VALUE) {
            cachedValue = MAX_VALUE;
        }

    }


    @Override
    public void prepareToNormalize(BaseInformationRecords.BaseInformation record, int indexOfRecord) {
        super.prepareToNormalize(record, indexOfRecord);
        cachedValue = 0;
        int index = 0;
        final List<BaseInformationRecords.CountInfo> countsList = sortedCountRecord.getSamples(this.sampleIndex).getCountsList();
        for (BaseInformationRecords.CountInfo count : countsList) {
            cachedValue |= (count.getIsCalled() ? 1 : 0) << index;
            index++;
            if (cachedValue > MAX_VALUE) {
                break;
            }
        }
        if (cachedValue > MAX_VALUE) {
            cachedValue = MAX_VALUE;
        }
        if (cachedValue == 0) {
            unableToRepresent.warn(LOG, "Unable to represent genotype, reached past max index (ploidy + extra-genotype)=" + index);
        }
    }

    @Override
    public float produceLabel(BaseInformationRecords.BaseInformation record, int labelIndex) {
        int n = numberOfLabels();
        float v = epsilon / (n - 1f);
        return (cachedValue == labelIndex) ? 1f - epsilon : v;
    }


    public static final String PLOIDY_PROPERTY = "genotypes.ploidy";

    @Override
    public void configure(Properties readerProperties) {

        String value = readerProperties.getProperty(PLOIDY_PROPERTY);
        try {
            maxCalledAlleles = Integer.parseInt(value) + 1;
        } catch (NumberFormatException e) {
            throw new RuntimeException("Unable to read ploidy from sbi properties file.");
        }
    }
}
