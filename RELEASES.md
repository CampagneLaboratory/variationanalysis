# SNAPSHOTS

Snapshot releases can be created as follows:

1. Make sure the correct version of Goby3 is built and installed in the local maven repository (typically do mvn install in the goby3 project.)
2. ````build-cpu.sh````
2. ```cd formal-releases```
3. ````./prepare-release.sh VERSION````
Where version is for instance 1.2-SNAPSHOT

# RELEASES

The same procedure is followed, but using a formal release number in the last step.

```
./prepare-release.sh 1.2
```