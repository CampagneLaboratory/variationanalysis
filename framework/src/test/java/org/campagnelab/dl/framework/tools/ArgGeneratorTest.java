package org.campagnelab.dl.framework.tools;

import org.junit.Test;

import static org.junit.Assert.*;

/**
 * Test the argument generator.
 * Created by fac2003 on 11/29/17.
 */
public class ArgGeneratorTest {
    @Test
    public void configureWithFlag() throws Exception {
        ArgGenerator generator = new ArgGenerator(1212);
        generator.configureWithString("--constant-learning-rates\n" +
                "flag\n" +
                "\n");
        assertEquals(
                "--constant-learning-rates\n" +
                        "\n" +
                        "--constant-learning-rates\n" +
                        "--constant-learning-rates\n" +
                        "--constant-learning-rates\n" +
                        "\n" +
                        "\n" +
                        "--constant-learning-rates\n" +
                        "--constant-learning-rates\n",
                generator.generateCommands(10));
    }

    @Test
    public void configureWithInt() throws Exception {
        ArgGenerator generator = new ArgGenerator(1212);
        generator.configureWithString("--num\n" +
                "int\n" +
                "1\n" +
                "10\n" +
                "\n");
        assertEquals(
                "--num 9\n" +
                        "--num 10\n" +
                        "--num 3\n" +
                        "--num 2\n" +
                        "--num 3\n" +
                        "--num 3\n" +
                        "--num 6\n" +
                        "--num 7\n" +
                        "--num 8\n" +
                        "--num 2",
                generator.generateCommands(10));
    }


    @Test
    public void configureMultiple() throws Exception {
        ArgGenerator generator = new ArgGenerator(1212);
        generator.configureWithString(
                "--num\n" +
                        "int\n" +
                        "1\n" +
                        "10\n" +
                        "\n" +
                        "--flag\n" +
                        "flag\n" +
                        "\n" +
                        "-x\n" +
                        "int\n" +
                        "100\n" +
                        "200\n" +
                        "" +
                        "\n");
        assertEquals(
                "--num 9  -x 135\n" +
                        "--num 2 --flag -x 184\n" +
                        "--num 6 --flag -x 122\n" +
                        "--num 2  -x 166\n" +
                        "--num 10 --flag -x 129\n" +
                        "--num 9 --flag -x 136\n" +
                        "--num 2 --flag -x 140\n" +
                        "--num 2  -x 116\n" +
                        "--num 10  -x 114\n" +
                        "--num 10 --flag -x 133",
                generator.generateCommands(10));
    }


}