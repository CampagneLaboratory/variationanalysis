package org.campagnelab.dl.genotype.mappers;

import org.campagnelab.dl.framework.mappers.ConfigurableFeatureMapper;
import org.campagnelab.dl.framework.mappers.OneHotBaseFeatureMapper;
import org.campagnelab.dl.framework.mappers.OneHotBaseLabelMapper;
import org.campagnelab.dl.genotype.helpers.GenotypeHelper;
import org.campagnelab.dl.varanalysis.protobuf.BaseInformationRecords;
import org.datavec.api.transform.transform.categorical.CategoricalToOneHotTransform;

import java.util.Properties;

/**
 * Label that encodes  the combination of called alleles in a sample using one-hot encoding and a softmax for decoding.<p/>
 * <p>
 * Created by fac2003 on 2/7/17.
 */
public class SoftmaxLabelMapper extends CountSortingLabelMapper implements ConfigurableFeatureMapper {


    private final float epsilon;
    public int ploidy;

    /**
     * @param sortCounts
     * @param ploidy
     * @param epsilon    amount of label smoothing to apply.
     */
    public SoftmaxLabelMapper(boolean sortCounts, int ploidy, float epsilon) {

        super(sortCounts);
        assert sortCounts : "SoftmaxLabelMapper requires sorted counts.";
        this.ploidy = ploidy;
        this.epsilon = epsilon;
    }

    public SoftmaxLabelMapper(boolean sortCounts, int ploidy) {

        this(sortCounts, ploidy, 0);
    }

    @Override
    public int numberOfLabels() {
        return (int) Math.pow(2, ploidy);
    }

    private int cachedValue;

    protected void setCachedValue(boolean... isCalled) {
        cachedValue = 0;
        int index = 0;
        for (boolean called : isCalled) {
            cachedValue |= (called ? 1 : 0) << index;
            index++;
        }
    }

    @Override
    public void prepareToNormalize(BaseInformationRecords.BaseInformation record, int indexOfRecord) {
        super.prepareToNormalize(record,indexOfRecord);
        cachedValue = 0;
        int index = 0;
        for (BaseInformationRecords.CountInfo count : sortedCountRecord.getSamples(0).getCountsList()) {
            cachedValue |= (count.getIsCalled() ? 1 : 0) << index;
            index++;
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
            ploidy = Integer.parseInt(value);
        } catch (NumberFormatException e) {
            throw new RuntimeException("Unable to read ploidy from sbi properties file.");
        }
    }
}
