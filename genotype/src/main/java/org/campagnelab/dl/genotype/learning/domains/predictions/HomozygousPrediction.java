package org.campagnelab.dl.genotype.learning.domains.predictions;

import org.campagnelab.dl.framework.domains.prediction.BinaryClassPrediction;
import org.campagnelab.dl.framework.domains.prediction.MultiClassPrediction;
import org.campagnelab.dl.genotype.mappers.HomozygousLabelsMapper;
import org.campagnelab.dl.varanalysis.protobuf.BaseInformationRecords;

/**
 * Created by rct66 on 11/12/16.
 */
public class HomozygousPrediction extends MultiClassPrediction {



    public <BaseInformation> void inspectRecord(BaseInformationRecords.BaseInformation currentRecord) {

        HomozygousLabelsMapper mapper = new HomozygousLabelsMapper();

        for (int i = 0; i < 11; i++){
            if (mapper.produceLabel(currentRecord,i) > 0.99F){
                trueClassIndex = i;
                return;
            }

        }
        trueClassIndex = -1;
        return;



    }
}

