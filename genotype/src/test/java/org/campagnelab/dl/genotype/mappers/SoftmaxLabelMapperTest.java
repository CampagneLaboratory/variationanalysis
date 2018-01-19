package org.campagnelab.dl.genotype.mappers;

import com.google.common.collect.Lists;
import it.unimi.dsi.fastutil.floats.FloatArrayList;
import org.junit.Test;

import static org.junit.Assert.assertEquals;

/**
 * Created by fac2003 on 12/20/16.
 */
public class SoftmaxLabelMapperTest {
    SoftmaxLabelMapper mapper = new SoftmaxLabelMapper(0,
            /*sort counts */false,
            /*ploidy */ 2,
            /*epsilon */ 0);

    @Test
    public void testDistinctAlleles() {
        mapper.maxCalledAlleles = 2;

        assertEquals("[1.0, 0.0, 0.0, 0.0, 0.0]", map(false, false).toString()); // 0 -> 0
        assertEquals("[0.0, 1.0, 0.0, 0.0, 0.0]", map(true, false).toString()); // 1 -> 1
        assertEquals("[0.0, 0.0, 1.0, 0.0, 0.0]", map(false, true).toString()); // 10 -> 2
        assertEquals("[0.0, 0.0, 0.0, 1.0, 0.0]", map(true, true).toString()); // 11 ->3

        mapper.maxCalledAlleles = 3;
        assertEquals("[1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]", map(false, false,false).toString()); // 000 -> 0
        assertEquals("[0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]", map(false, true, false).toString()); // 010 -> 2
        assertEquals("[0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0]", map(false, false, true).toString()); // 100 -> 4
        assertEquals("[0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0]", map(true, true, true).toString()); // 111 ->7

    }

    private FloatArrayList map(boolean... isCalled) {
        mapper.setCachedValue(isCalled);
        FloatArrayList result = new FloatArrayList();

        for (int featureLabel = 0; featureLabel < mapper.numberOfLabels(); featureLabel++) {
            result.add(mapper.produceLabel(null, featureLabel));
        }
        return result;
    }
}