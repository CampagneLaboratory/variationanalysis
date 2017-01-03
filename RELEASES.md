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
 2. Tag the release in git as "rVERSION"
 3. Update the tag's value in the scm section of the root pom.
 4. Deploy on the staging repository with ````mvn deploy -U````  (you need the CampagneLaboratory PGP keys for this step)
 5. Log in with the CampagneLaboratory at [Maven staging repository](https://oss.sonatype.org/#stagingRepositories) 
 6. Delete all the *bin-native* files from somatic and genotype modules in the CampagneLab repo
 7. Release the repository (select the repo and click on the release button above the list of repositories)