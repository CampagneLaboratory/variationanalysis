package org.campagnelab.dl.varanalysis.tools;

import org.apache.commons.io.FileUtils;
import org.campagnelab.dl.varanalysis.learning.TrainingArguments;
import org.junit.Test;

import java.io.File;
import java.io.IOException;

import static org.junit.Assert.assertEquals;

/**
 * Created by fac2003 on 11/5/16.
 */
public class ConditionRecordingToolTest {
    @Test
    public void tryTags() throws IOException {
        TrainingArguments args = new TrainingArguments();
        args.architectureClassname = "ABC";
        args.modelConditionFilename = "test-results/model-conditions/1.txt";
        args.learningRate = 1.0d;
        ConditionRecordingTool<TrainingArguments> tool = new ConditionRecordingTool<TrainingArguments>() {
            @Override
            public TrainingArguments createArguments() {
                return args;
            }

            @Override
            public TrainingArguments args() {
                return args;
            }

            @Override
            public void execute() {

            }
        };
        String commandLine = "--model-conditions test-results/model-conditions/1.txt -r 1.0 -t T -v V --random-seed 1478359791323 ";
        new File(args.modelConditionFilename).delete();
        tool.parseArguments(commandLine.split(" "), "name", args);
        tool.writeModelingConditions(args);
        assertEquals("Model condition file does not match expected.",
                FileUtils.readFileToString(new File(args.modelConditionFilename), "utf-8"),
                "Tag|Results|Specified_Arguments|Default_Arguments|Classname\n" +
                        "LH4CV8||--learning-rate 1.0 --model-conditions test-results/model-conditions/1.txt --random-seed 1478359791323 --training-sets [T] --validation-set V|--dropout-rate null --early-stopping-num-epochs 10 --error-enrichment false --experimental-condition not_specified --feature-mapper org.campagnelab.dl.model.utils.mappers.FeatureMapperV19 --learning-rate 1.0 --max-epochs 2147483647 --mini-batch-size 32 --model-conditions test-results/model-conditions/1.txt --net-architecture ABC --num-errors-added 16 --num-training 2147483647 --num-validation 2147483647 --parameter-precision FP32 --previous-model-name bestAUC --previous-model-path null --random-seed 1478359791323 --regularization-rate null --training-sets [T] --validate-every 1 --validation-set V|\n");
    }
}