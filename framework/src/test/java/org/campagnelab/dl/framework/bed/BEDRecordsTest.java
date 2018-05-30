package org.campagnelab.dl.framework.bed;

import org.junit.Test;

import java.io.StringReader;

import static org.junit.Assert.*;

public class BEDRecordsTest {
    @Test
    public void overlaps() throws Exception {

        String text_bed="1\t840550\t842913\n" +
                "1\t842947\t843250\n" +
                "1\t845279\t845442\n" +
                "1\t845563\t845638\n";
        BEDRecords records = BedLoader.loadBedFile(new StringReader(text_bed));
        assertTrue(records.overlaps("1",840550,840550));
        assertTrue(records.overlaps("1",840550,840551));
        assertTrue(records.overlaps("1",840550,842913));
        assertFalse(records.overlaps("1",842913,842913));

        //between first and second region:
        assertFalse(records.overlaps("1",842914,842946));

        // second region:
        assertTrue(records.overlaps("1",842947,842947));
        assertTrue(records.overlaps("1",842947,840551));
        assertTrue(records.overlaps("1",840550,843250));
        assertTrue(records.overlaps("1",840549,843251));
        assertFalse(records.overlaps("1",843251,843251));

    }

}