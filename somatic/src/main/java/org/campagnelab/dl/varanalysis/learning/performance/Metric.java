package org.campagnelab.dl.varanalysis.learning.performance;

public class Metric {
        String name;
        boolean largerIsBetter;

        public Metric(String name, boolean largerIsBetter) {
            this.name = name;
            this.largerIsBetter = largerIsBetter;
        }
    }