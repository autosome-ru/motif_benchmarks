#!/usr/bin/env sh
exec java $JAVA_OPTIONS -cp .:jstacs-2.3.jar PWMBench LOG "$@"
