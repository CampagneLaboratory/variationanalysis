package org.campagnelab.dl.framework.performance;

public class Metric {
        String name;
        boolean largerIsBetter;

        public Metric(String name, boolean largerIsBetter) {
            this.name = name;
            this.largerIsBetter = largerIsBetter;
        }
    }