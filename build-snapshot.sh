#!/usr/bin/env bash
echo "This script rebuilds the project using the latest DL4J snapshot"
mvn clean
mvn install -P CPU $*
mvn install -P GPU-SNAPSHOT $*