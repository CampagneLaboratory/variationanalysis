This module provides CUDA specific code to configure the project
when using CUDA PGUs. Code in this module is only included in packaged
JAR files when the GPU maven profile is specified when building (i.e.,
mvn -PGPU).

The code included in this module is executed when the program calls
```java
org.campagnelab.dl.framework.gpu.InitializeGpu.initialize()
```