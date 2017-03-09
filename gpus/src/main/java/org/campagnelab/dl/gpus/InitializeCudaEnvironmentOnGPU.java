package org.campagnelab.dl.gpus;
import org.nd4j.jita.conf.CudaEnvironment;
/**
 * Created by fac2003 on 12/1/16.
 */
public class InitializeCudaEnvironmentOnGPU {
    public void InitializeCudaEnvironmentOnGPU() {
        // temp workaround for backend initialization
        Nd4j.create(1);
        final long GB = 1024 * 1024 * 1024L;
        CudaEnvironment.getInstance().getConfiguration()
                .enableDebug(false)
                .allowMultiGPU(true)
                .setMaximumGridSize(1024)
                .setMaximumBlockSize(1024)
                .setMaximumDeviceCacheableLength(1*GB)
                .setMaximumDeviceCache(4L * GB)
                .setMaximumHostCacheableLength(1*GB)
                .setMaximumHostCache(16L *GB)
                // cross-device access is used for faster model averaging over pcie
                .allowCrossDeviceAccess(true);
        System.out.println("Configured CUDA environment.");
    }
}
