package org.campagnelab.dl.genotype.helpers;

import it.unimi.dsi.fastutil.objects.ObjectArraySet;
import org.campagnelab.goby.util.Variant;
import org.junit.Test;

import java.util.Set;

import static junit.framework.TestCase.assertFalse;
import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

/**
 * Created by fac2003 on 12/25/16.
 */
public class GenotypeHelperTest {
    @Test
    public void isNoCall() throws Exception {

    }

    @Test
    public void genotypeHasAlleleOrIndel() throws Exception {
        Set<Variant.FromTo> trueAlleles = new ObjectArraySet<>();
        trueAlleles.add(new Variant.FromTo("C", "C"));
        trueAlleles.add(new Variant.FromTo("CGGG", "C-GG"));
        assertTrue(GenotypeHelper.genotypeHasAlleleOrIndel(trueAlleles, "CGGGAA", "C-GGAA"));
        assertFalse(GenotypeHelper.genotypeHasAlleleOrIndel(trueAlleles, "CGGGA", "C-GGAA"));
        assertTrue(GenotypeHelper.genotypeHasAlleleOrIndel(trueAlleles, "C", "C"));
    }

    @Test
    public void testGenotypes() {
        assertEquals(true, GenotypeHelper.isVariant(true, "CTG/CTG", "C--"));
        assertEquals(false, GenotypeHelper.isVariant(false, "CTG/CTG", "C--"));
        assertEquals(true, GenotypeHelper.isVariant(true, "A--|A--", "ACA"));
        assertEquals(false, GenotypeHelper.isVariant(false, "A--|A--", "ACA"));
        assertEquals(true, GenotypeHelper.isNoCall("N|N"));
        assertEquals(true, GenotypeHelper.isNoCall("N/N"));
        assertEquals(false, GenotypeHelper.isNoCall("C|N"));
        assertEquals(true, GenotypeHelper.isVariant("CT|CT", "C-"));

        assertEquals(false, GenotypeHelper.isIndel("T", "A|C"));
        assertEquals(true, GenotypeHelper.isIndel("ACA", "A--|ACA"));
        assertEquals(true, GenotypeHelper.isIndel("G--", "G|T"));
        assertEquals(false, GenotypeHelper.isVariant("A", "A"));
        assertEquals(true, GenotypeHelper.isVariant("A|C", "C"));
        assertEquals(true, GenotypeHelper.isVariant("CA|C", "C"));
        assertEquals(true, GenotypeHelper.isVariant("CA|CT", "CT"));
        assertEquals(true, GenotypeHelper.isVariant("CA-|CT", "CT"));
        assertEquals(true, GenotypeHelper.isVariant(true, "C-|CT", "CT"));
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

    @Test
    public void checkPad() {
        assertEquals("ACA-", GenotypeHelper.pad("ACA", 4));
        assertEquals("GATA|G---", GenotypeHelper.padMulti("GATA|G", 4));


    }
}