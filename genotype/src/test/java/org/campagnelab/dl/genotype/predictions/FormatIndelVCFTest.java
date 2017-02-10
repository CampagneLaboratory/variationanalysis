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

    //A to: A/T  -> from: A to: A/T
    @Test
    public void formatIndelVCF3() throws Exception {
        String from = "A";
        Set<String> to = new ObjectArraySet<>();
        to.add("A");
        to.add("T");
        FormatIndelVCF format = new FormatIndelVCF(from,to, 'A');
        assertEquals("A",format.fromVCF);
        assertEquals(2,format.toVCF.size());
        assertTrue(format.toVCF.contains("A"));
        assertTrue(format.toVCF.contains("T"));
    }

    //A to: A/A  -> from: A to: A/A
    @Test
    public void formatIndelVCF4() throws Exception {
        String from = "A";
        Set<String> to = new ObjectArraySet<>();
        to.add("A");
        FormatIndelVCF format = new FormatIndelVCF(from,to, 'A');
        assertEquals("A",format.fromVCF);
        assertEquals(1,format.toVCF.size());
        assertTrue(format.toVCF.contains("A"));
    }

    //handle snp and del at same pos
    //from: G-C to: T/GTC  -> from: G to: T/GT
    @Test
    public void formatIndelVCF5() throws Exception {
        String from = "G-C";
        Set<String> to = new ObjectArraySet<>();
        to.add("T");
        to.add("GTC");
        FormatIndelVCF format = new FormatIndelVCF(from,to, 'G');
        assertEquals("G",format.fromVCF);
        assertEquals(2,format.toVCF.size());
        assertTrue(format.toVCF.contains("T"));
        assertTrue(format.toVCF.contains("GT"));

    }


    //handle snp and insertion at same pos
    //from: GAC to: T/G-C  -> from: GA to: TA/G
    @Test
    public void formatIndelVCF6() throws Exception {
        String from = "GAC";
        Set<String> to = new ObjectArraySet<>();
        to.add("T");
        to.add("G-C");
        FormatIndelVCF format = new FormatIndelVCF(from,to, 'G');
        assertEquals("GA",format.fromVCF);
        assertEquals(2,format.toVCF.size());
        assertTrue(format.toVCF.contains("TA"));
        assertTrue(format.toVCF.contains("G"));

    }



    //handle simple snp! don't want to accidentally change something that doesn't need fixing
    //from: A to: C/A  -> from: A  to: C/A
    @Test
    public void formatIndelVCF7() throws Exception {
        String from = "A";
        Set<String> to = new ObjectArraySet<>();
        to.add("A");
        to.add("C");
        FormatIndelVCF format = new FormatIndelVCF(from,to, 'A');
        assertEquals("A",format.fromVCF);
        assertEquals(2,format.toVCF.size());
        assertTrue(format.toVCF.contains("A"));
        assertTrue(format.toVCF.contains("C"));

    }
}