#!/usr/bin/env bash
mvn clean
mvn install -P CPU
mvn install -P GPU