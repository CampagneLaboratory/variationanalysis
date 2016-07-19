package org.campagnelab.dl.model.utils;

import it.unimi.dsi.fastutil.floats.FloatAVLTreeSet;
import it.unimi.dsi.fastutil.io.BinIO;

import java.io.File;
import java.io.IOException;

/**
 * Created by rct66 on 7/19/16.
 * Defines the class which calibrates model probability using some formula, such as bayes theorem.
 */
public abstract class CalcCalibrator {
    FloatAVLTreeSet plantedMutSet = new FloatAVLTreeSet();
    FloatAVLTreeSet allRecSet = new FloatAVLTreeSet();
    int totalExamples;
    String modelPath;

    public void observe(float modelProb, boolean isMut){
        totalExamples++;
        if (isMut) {
            plantedMutSet.add(modelProb);
        } else {
            allRecSet.add(modelProb);
        }
    }

    //call to save stats to disk
    public void save() throws IOException {
        File mutFile =  new File(modelPath + "/mutSet");
        File allFile = new File(modelPath + "/allSet");
        mutFile.createNewFile();
        allFile.createNewFile();
        BinIO.storeObject(plantedMutSet,mutFile);
        BinIO.storeObject(allRecSet,allFile);
    }

    //call to load stats from disk
    public void load() throws IOException, ClassNotFoundException {
        File mutFile =  new File(modelPath + "/mutSet");
        File allFile = new File(modelPath + "/allSet");
        plantedMutSet = (FloatAVLTreeSet)BinIO.loadObject(mutFile);
        allRecSet = (FloatAVLTreeSet)BinIO.loadObject(allFile);
    }

    public CalcCalibrator(String modelPath, boolean loadStats) throws IOException, ClassNotFoundException {
        this.modelPath = modelPath;
        if (loadStats){
            load();
        }
    }

    public abstract float calibrateProb(float modelProb);



}
