package org.campagnelab.dl.genotype.tools;

import it.unimi.dsi.fastutil.objects.ObjectAVLTreeSet;
import it.unimi.dsi.fastutil.objects.ObjectArraySet;
import org.junit.Test;

import java.util.Set;
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
        Set<String> to = new ObjectArraySet<>();
        to.add("G");
        to.add("C");
        //to = G/C
        assertEquals("0/1", PredictG.codeGT(to, "C", alts));
        to.remove("C");
        //to = G/G
        assertEquals("1", PredictG.codeGT(to, "C", alts));
        to.remove("G");
        to.add("C");
        //to = C/C
        assertEquals("0", PredictG.codeGT(to, "C", alts));
    }

}