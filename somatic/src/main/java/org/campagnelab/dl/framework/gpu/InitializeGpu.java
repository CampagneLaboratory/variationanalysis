package org.campagnelab.dl.framework.gpu;

/**
 * Created by fac2003 on 12/1/16.
 */
public class InitializeGpu {

    public static void initialize() {
        try {
            // this class is provided by the gpus maven module if we are building with CUDA.
            Class<?> initializeCudaEnvironment = Class.forName("org.campagnelab.dl.gpus.InitializeCudaEnvironment");
            if (initializeCudaEnvironment != null) {

                // calling the constructor initializes the environment.
                initializeCudaEnvironment.newInstance();
            }
        } catch (Exception e) {
            throw new RuntimeException("Unable to initialie GPU/CUDA environment.", e);
        }
    }
}
