package org.campagnelab.dl.varanalysis.learning.features;

import it.unimi.dsi.fastutil.floats.FloatArrayList;

/**
 * Intermediate representation of neural net features.
 * Created by fac2003 on 5/25/16.
 *
 * @author Fabien Campagne
 */
public class Features {
    FloatArrayList values;

    public Features(FloatArrayList values) {
        this.values = values;
    }

    public Features(Features another) {
        this.values = new FloatArrayList(another.values);
    }

    public Features(int numFeatures) {
        this.values = new FloatArrayList(numFeatures);
    }

    public float getFeatureValue(int index) {
        return values.getFloat(index);
    }
}
