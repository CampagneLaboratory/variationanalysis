package org.campagnelab.dl.genotype.helpers;

import org.junit.Test;

import static org.junit.Assert.assertEquals;

/**
 * Created by fac2003 on 12/25/16.
 */
public class GenotypeHelperTest {
    @Test
    public void testGenotypes() {
        assertEquals(false, GenotypeHelper.isVariant("CTG/CTG", "C"));
        assertEquals(false, GenotypeHelper.isVariant("CTG|CTG", "C"));
        assertEquals(true, GenotypeHelper.isNoCall("N|N"));
        assertEquals(true, GenotypeHelper.isNoCall("N/N"));
        assertEquals(false, GenotypeHelper.isNoCall("C|N"));
        assertEquals(false, GenotypeHelper.isVariant("CT|CT", "C"));

        assertEquals(false, GenotypeHelper.isIndel("A|C"));
        assertEquals(true, GenotypeHelper.isIndel("-|C"));
        assertEquals(true, GenotypeHelper.isIndel("-|T"));
        assertEquals(false, GenotypeHelper.isVariant("A", "A"));
        assertEquals(true, GenotypeHelper.isVariant("A|C", "C"));
        assertEquals(true, GenotypeHelper.isVariant("CA|C", "C"));
        assertEquals(true, GenotypeHelper.isVariant("CA|CT", "CT"));
        assertEquals(true, GenotypeHelper.isVariant("CA-|CT", "CT"));
        assertEquals(true, GenotypeHelper.isVariant(true,"C-|CT", "CT"));
        assertEquals(false, GenotypeHelper.isVariant(false, "C-|CT", "CT"));

        assertEquals(false, GenotypeHelper.isVariant("A", "A"));
        assertEquals(true, GenotypeHelper.isVariant("A", "C"));

    }
}