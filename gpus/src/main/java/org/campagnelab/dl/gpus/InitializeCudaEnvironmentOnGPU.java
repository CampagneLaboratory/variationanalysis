package org.campagnelab.dl.gpus;
import org.nd4j.jita.conf.CudaEnvironment;
/**
 * Created by fac2003 on 12/1/16.
 */
public class InitializeCudaEnvironmentOnGPU {
    public void InitializeCudaEnvironment() {
        CudaEnvironment.getInstance().getConfiguration()
                .enableDebug(false)
                .allowMultiGPU(true)
                .setMaximumGridSize(512)
                .setMaximumBlockSize(512)
                .setMaximumDeviceCacheableLength(1024 * 1024 * 1024L)
                .setMaximumDeviceCache(8L * 1024 * 1024 * 1024L)
                .setMaximumHostCacheableLength(1024 * 1024 * 1024L)
                .setMaximumHostCache(8L * 1024 * 1024 * 1024L)
                // cross-device access is used for faster model averaging over pcie
                .allowCrossDeviceAccess(true);
        System.out.println("Configured CUDA environment.");
    }
}
