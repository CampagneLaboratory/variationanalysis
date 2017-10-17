package org.campagnelab.dl.genotype.segments;

import org.junit.Test;

import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.OutputStream;
import java.util.Properties;

import static org.junit.Assert.*;

/**
 * Created by mas2182 on 10/17/17.
 */
public class SegmentLabelMapperTest {

    SegmentLabelMapper labelMapper = new SegmentLabelMapper(3);

    @Test(expected = IllegalArgumentException.class)
    public void mapWithError() throws Exception {
        labelMapper.map("A/G/T/A");
    }

    @Test
    public void map() throws Exception {
        float[] label = labelMapper.map("A/G/T");
        boolean found = false;
        for (float p : label) {
           if (p == 1L)
               if (!found)
                found = true;
           else
               assertTrue("More than one position is marked with 1",false);
        }
    }

    @Test
    public void writeMap() throws Exception {
        Properties props = new Properties();
        labelMapper.writeMap(props);
        assertTrue("Invalid number of properties mapped", props.size() == 36);
    }

    @Test
    public void numberOfLabels() throws Exception {
        assertEquals("Invalid number of label", 35, labelMapper.numberOfLabels());
    }

}