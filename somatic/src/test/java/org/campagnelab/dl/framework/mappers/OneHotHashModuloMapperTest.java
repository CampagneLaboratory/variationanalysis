package org.campagnelab.dl.framework.mappers;

import org.junit.Test;

import static org.junit.Assert.assertEquals;

/**
 * Created by fac2003 on 11/25/16.
 */
public class OneHotHashModuloMapperTest {
    @Test
    public void testObjectOneHot() {

        Object o[] = {"a s", 1, "F", 'g'};
        OneHotHashModuloMapper<Object> mapper;
        int[] numElements = {2, 10};
        for (int numElement : numElements) {
            System.out.printf("With %d elements:", numElement);
            mapper = new OneHotHashModuloMapper<>(numElement, object -> object);
            for (Object a : o) {
                mapper.prepareToNormalize(a, 0);
                for (int featureIndex = 0; featureIndex < mapper.numberOfFeatures(); featureIndex++) {
                    float feature = mapper.produceFeature(a, featureIndex);
                    System.out.println("object: " + a + " featureIndex: " + featureIndex + " feature: " + feature);
                    if (numElement == 2 && "a s".equals(a) && featureIndex == 0) {
                        assertEquals(1, feature, 0.01);
                    }
                    if (numElement == 2 && "a s".equals(a) && featureIndex == 1) {
                        assertEquals(0, feature, 0.01);
                    }

                    if (numElement == 10 && a.equals('g')) {
                        if (featureIndex == 3) {
                            assertEquals(1, feature, 0.01);
                        } else {
                            assertEquals(0, feature, 0.01);
                        }
                    }

                }
            }
        }


    }
}