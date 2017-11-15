package org.campagnelab.dl.genotype.segments;

import org.junit.Test;

import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.OutputStream;
import java.util.Arrays;
import java.util.List;
import java.util.Properties;

import static org.junit.Assert.*;

/**
 * Created by mas2182 on 10/17/17.
 */
public class SegmentLabelMapperTest {

    SegmentLabelMapper labelMapper = new SegmentLabelMapper(3);

    @Test(expected = IllegalArgumentException.class)
    public void mapWithError() throws Exception {
        labelMapper.map("A/G/T/A", Arrays.asList(0, 1, 2, 0));
    }

    @Test
    public void map() throws Exception {
        float[] mapped = labelMapper.map("A/G/T", Arrays.asList(0, 1, 2));
        System.out.println(Arrays.toString(mapped));
        boolean found = false;
        for (float p : mapped) {
           if (p == 1L)
               if (!found)
                found = true;
           else
               assertTrue("More than one position is marked with 1",false);
        }
        assertTrue("Unable to map",found);
    }

    @Test
    public void writeMap() throws Exception {
        Properties props = new Properties();
        labelMapper.writeMap(props);
        assertEquals("Invalid number of properties mapped", 217, props.size());

    }

    @Test
    public void numberOfLabels() throws Exception {
        assertEquals("Invalid number of label", 216, labelMapper.numberOfLabels());
    }

    @Test
    public void testSortByIndices1() throws Exception {
        String testString = "A";
        List<Integer> testIndices = Arrays.asList(9);
        assertEquals("A", SegmentLabelMapper.sortByIndices(testString, testIndices));
    }

    @Test
    public void testSortByIndices2() throws Exception {
        String testString = "ABC";
        List<Integer> testIndices = Arrays.asList(0, 1, 2);
        assertEquals("ABC", SegmentLabelMapper.sortByIndices(testString, testIndices));
    }

    @Test
    public void testSortByIndices3() throws Exception {
        String testString = "ABC";
        List<Integer> testIndices = Arrays.asList(3, 9, 6);
        assertEquals("ACB", SegmentLabelMapper.sortByIndices(testString, testIndices));
    }
}