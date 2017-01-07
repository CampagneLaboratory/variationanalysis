package org.campagnelab.dl.genotype.tools;

import it.unimi.dsi.fastutil.objects.ObjectAVLTreeSet;
import org.junit.Test;

import java.util.SortedSet;

import static org.junit.Assert.*;

/**
 * Created by fac2003 on 1/6/17.
 */
public class PredictGTest {
    @Test
    public void codeGT() throws Exception {
        SortedSet<String> alts = new ObjectAVLTreeSet<>();
        alts.add("G");
        assertEquals("0/1", PredictG.codeGT("G/C", "C", alts));
        assertEquals("1", PredictG.codeGT("G/G", "C", alts));
        assertEquals("0", PredictG.codeGT("C/C", "C", alts));
    }

}