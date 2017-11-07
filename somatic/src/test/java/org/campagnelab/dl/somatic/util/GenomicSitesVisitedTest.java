package org.campagnelab.dl.somatic.util;

import org.junit.Test;

import static org.junit.Assert.*;

public class GenomicSitesVisitedTest {

    @Test
    public void testVisit() {
        GenomicSitesVisited vis = new GenomicSitesVisited();

        vis.visit(0, 1);
        assertTrue(vis.wasVisited(0, 1));
        assertFalse(vis.wasVisited(0, 2));
        assertFalse(vis.wasVisited(1, 0));
        vis.visit(0xFFFF, 0);
        assertFalse(vis.wasVisited(0, 0xFFFF));
        assertTrue(vis.wasVisited(0xFFFF, 0));
    }
}