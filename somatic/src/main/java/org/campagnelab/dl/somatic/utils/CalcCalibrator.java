package org.campagnelab.dl.somatic.utils;

import it.unimi.dsi.fastutil.floats.FloatAVLTreeSet;
import it.unimi.dsi.fastutil.io.BinIO;

import java.io.File;
import java.io.IOException;

/**
 * Created by rct66 on 7/19/16.
 * Defines the class which calibrates model probability using some formula, such as bayes theorem.
 *
 */
public abstract class CalcCalibrator {
    /**
     * These are two sorted sets containing the model probabilities encountered in PredictMutations (mutated and unmutated examples seperate).
     * They structured for easy access to the number of model probabilities greater than x, (or less than), for any given x.
     */
    FloatAVLTreeSet plantedMutSet = new FloatAVLTreeSet();
    FloatAVLTreeSet unMutSet = new FloatAVLTreeSet();
    int totalExamples;
    String modelPath;
    String prefix;

    public void observe(float modelProb, boolean isMut){
        totalExamples++;
        if (isMut) {
            plantedMutSet.add(modelProb);
        } else {
            unMutSet.add(modelProb);
        }
    }

    //call to save stats to disk
    public void save() throws IOException {
        File mutFile =  new File(modelPath + "/" + prefix + "mutSet");
        File unMutFile = new File(modelPath + "/" + prefix + "unMutSet");
        mutFile.createNewFile();
        unMutFile.createNewFile();
        BinIO.storeObject(plantedMutSet,mutFile);
        BinIO.storeObject(unMutSet,unMutFile);
    }

    //call to load stats from disk
    public void load() throws IOException, ClassNotFoundException {
        File mutFile =  new File(modelPath + "/" + prefix + "mutSet");
        File unMutFile = new File(modelPath + "/" + prefix + "unMutSet");
        plantedMutSet = (FloatAVLTreeSet)BinIO.loadObject(mutFile);
        unMutSet = (FloatAVLTreeSet)BinIO.loadObject(unMutFile);
    }

    public CalcCalibrator(String modelPath, String prefix, boolean loadStats) throws IOException, ClassNotFoundException {
        this.modelPath = modelPath;
        this.prefix = prefix;
        if (loadStats){
            load();
        }
    }

    public abstract float calibrateProb(float modelProb);



}
