package org.campagnelab.dl.genotype.mappers;

import org.campagnelab.dl.framework.mappers.FeatureCalculator;
import org.campagnelab.dl.somatic.genotypes.BaseGenotypeCountFactory;
import org.campagnelab.dl.somatic.genotypes.GenotypeCountFactory;
import org.campagnelab.dl.somatic.mappers.AbstractFeatureMapper;
import org.campagnelab.dl.somatic.mappers.GenotypeCount;
import org.campagnelab.dl.varanalysis.protobuf.BaseInformationRecords;
import org.nd4j.linalg.api.ndarray.INDArray;

/**
 * This is a simple feature mapper.
 * Each sample has the following information:
 * <pre> 69033640	11	false
 * position=14521   referenceIndex=0       isMutated=false
 * sample 0 counts=[10, 0, 0, 63, 0, 4, 0, 1, 62, 0]
 * <p>
 * Using these data, we can normalize counts and map them
 *
 * Created by rct66 on 11/16/16.
 *
 * @author Remi Torracinta
 */

public class SingleGenoTypeCountMapper extends AbstractFeatureMapperStripped<BaseInformationRecords.BaseInformationOrBuilder> {


    int sample;
    int genotype;
    boolean getForwardStrand;


    public SingleGenoTypeCountMapper(int sample, int genotype, boolean getForwardStrand) {
        this.sample = sample;
        this.genotype = genotype;
        this.getForwardStrand = getForwardStrand;
    }

    @Override
    public int numberOfFeatures() {
        return 1;
    }


    public void prepareToNormalize(BaseInformationRecords.BaseInformationOrBuilder record, int indexOfRecord) {
        //normalization will take place after concating this mapper with other count mappers
    }

    public float produceFeature(BaseInformationRecords.BaseInformationOrBuilder record, int featureIndex) {
        return produceFeatureInternal(record,featureIndex);
    }

    @Override
    public String getFeatureName(int featureIndex) {
        String sampleString = Integer.toString(sample)+"Sample";
        String genotypeString;
        switch (genotype) {
            case 0: genotypeString = "A"; break;
            case 1: genotypeString = "T"; break;
            case 2: genotypeString = "C"; break;
            case 3: genotypeString = "G"; break;
            case 4: genotypeString = "N"; break;
            case 5: genotypeString = "Indel1"; break;
            case 6: genotypeString = "Indel2"; break;
            case 7: genotypeString = "Indel3"; break;
            case 8: genotypeString = "Indel4"; break;
            default: genotypeString = "Indel5";
        }
        return ("sample"+sampleString+"Genotype"+genotypeString+(getForwardStrand?"Forward":"Backward"));
    }

    private float normalize(float value, int normalizationFactor) {
        //we are not normalizing here
        return value;
    }


    public float produceFeatureInternal(BaseInformationRecords.BaseInformationOrBuilder record, int featureIndex) {
        assert (featureIndex >= 0 && featureIndex < 1) : "This mapper only outputs 1 feature corresponding to one base count";
        BaseInformationRecords.CountInfo genoInfo;
        try {
            genoInfo = record.getSamples(sample).getCounts(genotype);
        } catch (NullPointerException|IndexOutOfBoundsException e){
            return 0;
        }
        float val = getForwardStrand?genoInfo.getGenotypeCountForwardStrand():genoInfo.getGenotypeCountReverseStrand();
        return val;
    }

}
