package org.campagnelab.dl.genotype.tools;


import org.campagnelab.dl.framework.mappers.ConfigurableFeatureMapper;
import org.campagnelab.dl.framework.mappers.FeatureNameMapper;
import org.campagnelab.dl.framework.tools.arguments.AbstractTool;

import java.io.IOException;
import java.io.StringReader;
import java.util.Properties;


/*
 *
 * @author rct66
 */
public class FeatureNamePrinter extends AbstractTool<FeatureNameArguments> {

    FeatureNameMapper mapper;

    public static void main(String[] args) {
        FeatureNamePrinter printer = new FeatureNamePrinter();
        printer.parseArguments(args, "FeatureNamePrinter", printer.createArguments());
        printer.execute();

    }

    public FeatureNamePrinter() {

    }


    @Override
    public FeatureNameArguments createArguments() {
        return new FeatureNameArguments();
    }

    @Override
    public void execute() {
        try {
            Class clazz = Class.forName(args().featureMapperClassname);
            mapper = (FeatureNameMapper) clazz.newInstance();
            Properties prop = new Properties();
            prop.load(new StringReader(
                    "stats.baseQuality.reverse.min=0.0\n" +
                    "numRecords=3224437\n" +
                    "stats.readMappingQuality.forward.min=0.0\n" +
                    "stats.baseQuality.forward.max=127.0\n" +
                    "stats.numVariationsInRead.max=170.0\n" +
                    "stats.genomicContextSize.min=21.0\n" +
                    "stats.readMappingQuality.reverse.min=0.0\n" +
                    "stats.baseQuality.reverse.max=127.0\n" +
                    "stats.insertSizes.min=-2.8001616E7\n" +
                    "stats.readMappingQuality.forward.max=70.0\n" +
                    "stats.genomicContextSize.max=21.0\n" +
                    "stats.readMappingQuality.reverse.max=70.0\n" +
                    "goby.version=\n" +
                    "stats.baseQuality.forward.min=0.0\n" +
                    "stats.insertSizes.max=2.8001616E7\n" +
                    "stats.numVariationsInRead.min=0.0"));
            ((ConfigurableFeatureMapper)mapper).configure(prop);
            for (int i = 0; i < mapper.numberOfFeatures(); i++){
                System.out.println(mapper.getFeatureName(i) + "\t");
            }
        } catch (IllegalAccessException e) {
            e.printStackTrace();
        } catch (InstantiationException e) {
            e.printStackTrace();
        } catch (ClassNotFoundException e) {
            e.printStackTrace();
        } catch (IOException e) {
            e.printStackTrace();
        }
        return;
    }
}