package org.campagnelab.dl.somatic.mappers;

import it.unimi.dsi.fastutil.doubles.DoubleArrayList;
import org.campagnelab.dl.varanalysis.protobuf.BaseInformationRecords;
import org.junit.Test;

import java.util.ArrayList;
import java.util.List;

import static org.junit.Assert.assertArrayEquals;
import static org.junit.Assert.assertEquals;

/**
 * Created by rct66 on 12/28/16.
 */
public class BamFlagMapperTest {

    @Test
    public void mapBits() throws Exception {
        final List<BaseInformationRecords.NumberWithFrequency> list = new ArrayList<>();
        BaseInformationRecords.NumberWithFrequency.Builder builder = BaseInformationRecords.NumberWithFrequency.newBuilder();
        builder.setNumber(1187);
        builder.setFrequency(1);
        list.add(builder.build());
        BamFlagMapper mapper = new BamFlagMapper(12, baseInformationOrBuilder -> list);
        //   mapper.prepareToNormalize(null, 0);
        double unit = 1;
        //     assertEquals(DoubleArrayList.wrap(new double[]{unit/5, unit/5, 0, 0, 0, unit/5, 0, unit/5, 0, 0, unit/5, 0}), DoubleArrayList.wrap(mapper.propFractions));

        builder = BaseInformationRecords.NumberWithFrequency.newBuilder();
        builder.setNumber(1187 + 4);
        builder.setFrequency(1);
        list.add(builder.build());
        mapper = new BamFlagMapper(12, baseInformationOrBuilder -> list);
        double d = 5 * 2 + 1;
        mapper.prepareToNormalize(null, 0);
        assertEquals(DoubleArrayList.wrap(new double[]{2 * unit / d, 2 * unit / d, unit / d, 0, 0, 2 * unit / d, 0, 2 * unit / d, 0, 0, 2 * unit / d, 0}), DoubleArrayList.wrap(mapper.propFractions));

    }
}