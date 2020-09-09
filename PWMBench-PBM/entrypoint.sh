#!/usr/bin/env sh
exec java $JAVA_OPTIONS -cp /app:/app/jstacs-2.3.jar:/app/json-simple-1.1.1.jar PWMBench "$@"
