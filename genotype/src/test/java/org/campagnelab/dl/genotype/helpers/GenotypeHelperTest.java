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
    @Test
    public void checkMatchingGenotypes() {
        assertEquals(true, GenotypeHelper.matchingGenotypesWithN("A", "N"));
        assertEquals(false, GenotypeHelper.matchingGenotypes("A", "N"));

        assertEquals(false, GenotypeHelper.matchingGenotypes("A", "C"));
        assertEquals(false, GenotypeHelper.matchingGenotypes("A", ""));
        assertEquals(true, GenotypeHelper.matchingGenotypesWithN("N", "A"));
        assertEquals(false, GenotypeHelper.matchingGenotypes("N", "A"));
        assertEquals(true, GenotypeHelper.matchingGenotypes("GT/GT", "G"));
        assertEquals(false, GenotypeHelper.matchingGenotypes("CCTA/CCTA", "A"));
        assertEquals(true, GenotypeHelper.matchingGenotypes("CCTA/CCTA", "C"));
        assertEquals(true, GenotypeHelper.matchingGenotypes("CCTA/CCTA", "CC"));
        assertEquals(true, GenotypeHelper.matchingGenotypes("CCTA/CCTA", "CCT"));
        assertEquals(false, GenotypeHelper.matchingGenotypes("CCTA/CCTA", "CCA"));
    }

    @Test
    public void checkHasGenotype() {
        assertEquals(false, GenotypeHelper.genotypeHasAllele("A/C", "N"));
        assertEquals(true, GenotypeHelper.genotypeHasAllele("A/C", "A"));
        assertEquals(true, GenotypeHelper.genotypeHasAllele("A/C", "C"));
        assertEquals(false, GenotypeHelper.genotypeHasAllele("A/C", "T"));
        assertEquals(false, GenotypeHelper.genotypeHasAllele("A/C", "T"));
        assertEquals(true, GenotypeHelper.genotypeHasAllele("GT/GT", "G"));


    }
}