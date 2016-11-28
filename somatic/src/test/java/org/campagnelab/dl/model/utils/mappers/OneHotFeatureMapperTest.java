package org.campagnelab.dl.model.utils.mappers;

import org.campagnelab.dl.varanalysis.protobuf.BaseInformationRecords;
import org.junit.Test;
import org.campagnelab.dl.somatic.mappers.BinaryFeatureMapper;

import static org.junit.Assert.assertEquals;

/**
 * Created by fac2003 on 7/12/16.
 */
public class OneHotFeatureMapperTest {
    @Test
    public void testOneHot() {
        BinaryFeatureMapper mapper;

        mapper=new BinaryFeatureMapper() {
            @Override
            public int getIntegerValue(BaseInformationRecords.BaseInformationOrBuilder record) {
                return 0;
            }
        };
        mapper.prepareToNormalize(null,0);
        for (int i=0;i<32;i++) {
            assertEquals(0f,mapper.produceFeature(null,i),0.1);
        }

        mapper=new BinaryFeatureMapper() {
            @Override
            public int getIntegerValue(BaseInformationRecords.BaseInformationOrBuilder record) {
                return 1;
            }
        };
        mapper.prepareToNormalize(null,0);
        for (int i=0;i<32;i++) {
            if (i!=0) {
                assertEquals(0f, mapper.produceFeature(null, i), 0.1);
            }
        }
        assertEquals(1f, mapper.produceFeature(null, 0), 0.1);

        mapper=new BinaryFeatureMapper() {
            @Override
            public int getIntegerValue(BaseInformationRecords.BaseInformationOrBuilder record) {
                return 0xFFFFFFFF;
            }
        };
        mapper.prepareToNormalize(null,0);
        for (int i=0;i<32;i++) {
            assertEquals(1f,mapper.produceFeature(null,i),0.1);
        }
    }

}