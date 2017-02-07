package org.campagnelab.dl.genotype.predictions;

import it.unimi.dsi.fastutil.objects.ObjectAVLTreeSet;
import it.unimi.dsi.fastutil.objects.ObjectArraySet;
import org.campagnelab.dl.genotype.tools.PredictG;
import org.junit.Test;

import java.util.Set;
import java.util.SortedSet;

import static org.junit.Assert.*;

/**
 * Created by rct66 on 2/7/17.
 */
public class FormatIndelVCFTest {

    //IE: from: GTAC to: G--C/G-AC -> from: GTA to: G/GA
    @Test
    public void formatIndelVCF() throws Exception {


        String from = "GTAC";
        Set<String> to = new ObjectArraySet<>();
        to.add("G--C");
        to.add("G-AC");
        FormatIndelVCF format = new FormatIndelVCF(from,to,'G');
        assertEquals("GTA",format.fromVCF);
        assertEquals(2,format.toVCF.size());
        assertTrue(format.toVCF.contains("G"));
        assertTrue(format.toVCF.contains("GA"));
    }

    //TGG to: T/T-G  -> from: TG to: TG/T
    @Test
    public void formatIndelVCF2() throws Exception {
        String from = "TGG";
        Set<String> to = new ObjectArraySet<>();
        to.add("T");
        to.add("T-G");
        FormatIndelVCF format = new FormatIndelVCF(from,to, 'T');
        assertEquals("TG",format.fromVCF);
        assertEquals(2,format.toVCF.size());
        assertTrue(format.toVCF.contains("TG"));
        assertTrue(format.toVCF.contains("T"));
    }
}