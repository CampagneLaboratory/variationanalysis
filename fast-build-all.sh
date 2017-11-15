#!/usr/bin/env bash
mvn install -P CPU -DskipTests=true $*
mvn package -P GPU -DskipTests=true $*