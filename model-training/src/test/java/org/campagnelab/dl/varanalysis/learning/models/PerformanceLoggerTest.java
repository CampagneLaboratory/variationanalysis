package org.campagnelab.dl.varanalysis.learning.models;

import org.junit.Test;

/**
 * Created by fac2003 on 7/24/16.
 */
public class PerformanceLoggerTest {
    PerformanceLogger logger = new PerformanceLogger("test-results/models/1213121212");

    @Test
    public void logWrite() throws Exception {
        logger.log("best", 0, 0, Float.NaN, Float.NaN);
        logger.log("best", 1, 0, Float.NaN, Float.NaN);
        logger.log("best", 1000, 1, 0.3f, Float.NaN);
        logger.log("best", 10000, 10, 0.3f, 0.5f);
        logger.write("best");
    }

    @Test
    public void logWriteSeveral() throws Exception {
        logger.log("best", 0, 0, Float.NaN, Float.NaN);
        logger.log("bestAUC", 1000, 1, 0.3f, 0.9f);
        logger.log("final", 10000, 10, 0.3f, 0.5f);
        logger.write();
    }
    @Test
    public void logWriteConditionId() throws Exception {
        PerformanceLogger    logger = new PerformanceLogger("test-results/models/conditionId");
        logger.setCondition("Condition1");
        logger.log("best", 0, 0, Float.NaN, Float.NaN);
        logger.log("best", 1, 0, Float.NaN, Float.NaN);
        logger.log("best", 1000, 1, 0.3f, Float.NaN);
        logger.log("best", 10000, 10, 0.3f, 0.5f);
        logger.write("best");
    }


}