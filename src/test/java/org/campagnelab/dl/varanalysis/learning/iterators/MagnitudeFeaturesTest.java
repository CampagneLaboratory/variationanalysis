package org.campagnelab.dl.varanalysis.learning.iterators;

import org.junit.Test;

import static org.junit.Assert.*;

/**
 * Created by fac2003 on 5/24/16.
 */
public class MagnitudeFeaturesTest {
    @Test
    public void mapFeatures() throws Exception {
        MagnitudeFeatures magnitude = new MagnitudeFeatures();

        assertEquals(31, magnitude.numberOfFeatures());
        // TODO complete the test once we have a PB class for record.

    }

}