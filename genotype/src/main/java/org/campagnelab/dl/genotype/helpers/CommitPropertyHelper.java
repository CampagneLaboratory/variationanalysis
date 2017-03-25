package org.campagnelab.dl.genotype.helpers;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.io.IOException;
import java.io.InputStream;
import java.util.Properties;

/**
 * Helper class to retrive commit info and store it to a Properties instance.
 * Created by fac2003 on 3/25/17.
 */
public class CommitPropertyHelper {
    private static final Logger LOG = LoggerFactory.getLogger(CommitPropertyHelper.class);

    public static void appendCommitInfo(Object somethingFromProject, String path, Properties destination) {

        InputStream in = somethingFromProject.getClass().getResourceAsStream(path);
        if (in == null) {
            LOG.error("Goby commit properties file ("+path+") not found in classpath. Unable to write info to sbip.");
        } else {
            Properties gobyProperties = new Properties();
            try {
                gobyProperties.load(in);
            } catch (IOException e) {
                LOG.error("Unable to load properties from  ("+path+").");
            }
            destination.putAll(gobyProperties);
        }

    }
}
