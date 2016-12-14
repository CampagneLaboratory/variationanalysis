## How to compile

To build the project [Maven](https://maven.apache.org) 3.2+ is required. Maven is a Java tool, so you must have [Java 1.8+](http://www.oracle.com/technetwork/java/javase/downloads/index.html) installed in order to proceed.

To check your local Maven installation, type the following in a terminal or in a command prompt:
```sh
mvn --version
```
If Maven is installed, it should print out your installed version and configuration:
```sh
Apache Maven 3.0.5 (r01de14724cdef164cd33c7c8c2fe155faf9602da; 2013-02-19 14:51:28+0100)
Maven home: /usr/local/Cellar/maven/3.0.5/libexec
Java version: 1.8.0_65, vendor: Oracle Corporation
Java home: /Library/Java/JavaVirtualMachines/jdk1.8.0_65.jdk/Contents/Home/jre
Default locale: en_US, platform encoding: UTF-8
OS name: "mac os x", version: "10.10.5", arch: "x86_64", family: "mac"
```
If the version is less than 3.2 or Maven is not installed, you need to install a new version of the tool.

### Upgrade your Maven installation

To install a new version of Maven, first, [download](https://maven.apache.org/download.html) the software and then follow the [installation instructions](https://maven.apache.org/install.html).

After the installation, check that the new version has been correctly installed:

```sh
mvn -version
Apache Maven 3.3.9 (bb52d8502b132ec0a5a3f4c09453c07478323dc5; 2015-11-10T11:41:47-05:00)
Maven home: /usr/local/Cellar/maven/3.3.9/libexec
Java version: 1.8.0_65, vendor: Oracle Corporation
Java home: /Library/Java/JavaVirtualMachines/jdk1.8.0_65.jdk/Contents/Home/jre
Default locale: en_US, platform encoding: UTF-8
OS name: "mac os x", version: "10.10.5", arch: "x86_64", family: "mac"
```

### Compile from the command line
In a terminal or in a command prompt, change to the project root folder (the folder where you cloned this git repo) and type the following command:
```sh
mvn package
```
The command line will print out various actions, and end with the following:
```sh
[INFO] ------------------------------------------------------------------------
[INFO] Reactor Summary:
[INFO] 
[INFO] Variation Analysis ................................. SUCCESS [  0.635 s]
[INFO] framework .......................................... SUCCESS [  4.780 s]
[INFO] somatic ............................................ SUCCESS [ 14.967 s]
[INFO] genotype ........................................... SUCCESS [  1.473 s]
[INFO] gpus ............................................... SUCCESS [  0.052 s]
[INFO] ------------------------------------------------------------------------
[INFO] BUILD SUCCESS
[INFO] ------------------------------------------------------------------------
[INFO] Total time: 22.600 s
[INFO] Finished at: 2016-12-09T16:43:41-05:00
[INFO] Final Memory: 363M/1223M
[INFO] ------------------------------------------------------------------------
```

You may skip unit tests with the command 

```sh
mvn clean package -DskipTests
```

### GPU support

You can compile the project for CUDA GPUs by activating the GPU profile:
```sh
mvn clean package -PGPU
```

This will create jar files that include CUDA support for DL4J and activate some GPU specific customizations.

### Compile in IntelliJ Idea

If building from within [IntelliJ Idea](https://www.jetbrains.com/idea/), make sure the IDE is configured to build with Maven 3.2+.

To check/change the Maven Integration settings:

- go to Preferences...>Build,Execution,Deployment>Build Tools and select Maven
- in the dialog, check the Version number reported under "Maven home directory"
- if the version is less than 3.2+, click on the browse button next to the input field and select the "Maven home" folder reported by the mvn --version command

### Troubleshooting

If you get an error similar to the following one during the build, it's because you are using an outdated Maven version:

```sh
[INFO] ------------------------------------------------------------------------
[INFO] Detecting the operating system and CPU architecture
[INFO] ------------------------------------------------------------------------
[INFO] os.detected.name: osx
[INFO] os.detected.arch: x86_64
[INFO] os.detected.version: 10.11
[INFO] os.detected.version.major: 10
[INFO] os.detected.version.minor: 11
[INFO] os.detected.classifier: osx-x86_64
[WARNING] Failed to inject repository session properties.
java.lang.NoSuchMethodError: org.apache.maven.execution.MavenSession.getRepositorySession()Lorg/eclipse/aether/RepositorySystemSession;
	at kr.motd.maven.os.RepositorySessionInjector.injectRepositorySession(RepositorySessionInjector.java:22)
```
In this case, review the installation according to what reported above and make sure that you are using Maven 3.2+.


