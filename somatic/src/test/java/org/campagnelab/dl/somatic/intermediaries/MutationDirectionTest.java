package org.campagnelab.dl.somatic.intermediaries;

import org.junit.Test;

import java.util.Arrays;

import static org.junit.Assert.assertEquals;

/**
 * Test the mutator on some specific examples.
 * Created by fac2003 on 5/27/16.
 */
public class MutationDirectionTest {
    @Test
    public void somaticFrequency() {
        double delta = 0.75;
        FirstSimulationStrategy strategy = new FirstSimulationStrategy(delta, delta, 0.1, 1);
        FirstSimulationStrategy.MutationDirection dir = null;

        dir = strategy.dirFromCounts(0, new int[]{50, 0, 0, 50});
        assertEquals("somaticFrequency must match",delta, dir.somaticFrequency, 0.1);

        dir = strategy.dirFromCounts(0,new int[]{80, 20, 0, 0});
        assertEquals("somaticFrequency must match",delta, dir.somaticFrequency,0.1 );

        dir = strategy.dirFromCounts(0, new int[]{80, 20, 3, 1});
        assertEquals("somaticFrequency must match",delta, dir.somaticFrequency,0.1);

        dir = strategy.dirFromCounts(0,new int[]{90, 1, 1, 1});
        assertEquals("somaticFrequency must match",delta , dir.somaticFrequency, 0.1);

        dir = strategy.dirFromCounts(0,new int[]{0, 3, 1, 100});
        assertEquals("somaticFrequency must match", delta , dir.somaticFrequency, 0.1);
    }

    @Test
    public void mutateTest() throws Exception {
        FirstSimulationStrategy strategy = new FirstSimulationStrategy(1, 1, 0.1, 1);
        final int mostCount = 100;
        final int secondMostCount = 12;
        //homozygous first:
        for (int numGenos = 5; numGenos < 8; numGenos++) {
            for (int maxIndex = 0; maxIndex < numGenos; maxIndex++) {
                if (maxIndex == 4) {
                    continue;
                }
                int[] counts = new int[numGenos];
                counts[maxIndex] = mostCount;
                for (int trial = 0; trial < 20; trial++) {
                    FirstSimulationStrategy.MutationDirection dir = strategy.dirFromCounts(0,counts);
                    System.out.println(Arrays.toString(counts) + " " + dir.oldBase + " -> " + dir.newBase);
                    assert (dir.newBase != dir.oldBase);
                    assert (dir.newBase != maxIndex);
                    assert (dir.newBase != 4);
                    assert (dir.newBase < numGenos);
                }
            }
        }
        System.out.println("\nHeterozygous:");
        //now heterozygous
        for (int numGenos = 5; numGenos < 8; numGenos++) {
            for (int maxIndex = 0; maxIndex < numGenos; maxIndex++) {
                if (maxIndex == 4) {
                    continue;
                }
                for (int secondMostIndex = 0; secondMostIndex < numGenos; secondMostIndex++) {
                    if (secondMostIndex == maxIndex || secondMostIndex == 4) {
                        continue;
                    }
                    int[] counts = new int[numGenos];
                    counts[maxIndex] = mostCount;
                    counts[secondMostIndex] = secondMostCount;
                    for (int trial = 0; trial < 20; trial++) {
                        FirstSimulationStrategy.MutationDirection dir = strategy.dirFromCounts(0,counts);
                        System.out.println(Arrays.toString(counts) + " " + dir.oldBase + " -> " + dir.newBase);
                        assert (dir.newBase != dir.oldBase);
                        assert (dir.newBase != maxIndex);
                        assert (dir.newBase != 4);
                        assert (dir.newBase != secondMostIndex);
                        assert (dir.newBase < numGenos);
                    }
                }

            }
        }

    }


}