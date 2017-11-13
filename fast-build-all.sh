#!/usr/bin/env bash
mvn mvn install -P CPU -DskipTests=true $*
mvn install -P GPU -DskipTests=true $*