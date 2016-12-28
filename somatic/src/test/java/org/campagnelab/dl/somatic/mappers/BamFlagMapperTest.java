package org.campagnelab.dl.somatic.mappers;

import org.junit.Test;

import static org.junit.Assert.*;

/**
 * Created by rct66 on 12/28/16.
 */
public class BamFlagMapperTest {
    @Test
    public void decodeProps() throws Exception {
        BamFlagMapper mapper = new BamFlagMapper(0);
        boolean[] decoded = mapper.decodeProps(1187);
        assertArrayEquals(new boolean[]{true,true,false,false,false,true,false,true,false,false,true,false},decoded);
    }

}