package org.campagnelab.dl.genotype.mappers;

import com.google.common.primitives.UnsignedLongs;
import it.unimi.dsi.fastutil.floats.FloatArrayList;
import org.junit.Test;

import static org.junit.Assert.*;

/**
 * Created by fac2003 on 12/20/16.
 */
public class NumDistinctAllelesLabelMapperTest {
    NumDistinctAllelesLabelMapper mapper=new NumDistinctAllelesLabelMapper(false,0);

    @Test
    public void testDistinctAlleles() {
        mapper.ploidy=2;
        assertEquals("[0.0, 1.0]",map("A/B").toString());
        assertEquals("[1.0, 0.0]",map("A/A").toString());
        assertEquals("[1.0, 0.0]",map("B/B").toString());

        mapper.ploidy=3;
        assertEquals("[0.0, 1.0, 0.0]",map("A/B").toString());
        assertEquals("[1.0, 0.0, 0.0]",map("A").toString());
        assertEquals("[0.0, 0.0, 1.0]",map("A/B/C").toString());

    }

    private FloatArrayList map(String trueGenotype) {
        FloatArrayList result=new FloatArrayList();
        for (int featureLabel=0;featureLabel<mapper.ploidy;featureLabel++) {
            result.add(mapper.label(featureLabel, trueGenotype));
        }
        return result;
    }
}