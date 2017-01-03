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

# MAVEN RELEASES
 1. Make sure the project version does not end with -SNAPSHOT in the poms
 2. Deploy on the staging repository with ````mvn deploy -U````  (you need the CampagneLaboratory PGP keys for this step)
 3. Log in with the CampagneLaboratory at [Maven staging repository](https://oss.sonatype.org/#stagingRepositories) 
 4. Delete all the *bin-native* files from somatic and genotype modules in the CampagneLab repo
 5. Release the repository ()