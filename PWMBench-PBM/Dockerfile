FROM openjdk:8-jdk-alpine
WORKDIR /app
COPY entrypoint.sh /app/entrypoint.sh
COPY PWMBench.java /app/PWMBench.java
COPY jstacs-2.3.jar /app/jstacs-2.3.jar
COPY json-simple-1.1.1.jar /app/json-simple-1.1.1.jar
RUN javac -cp .:jstacs-2.3.jar:json-simple-1.1.1.jar PWMBench.java

#########################

FROM openjdk:8-jre-alpine
LABEL maintainer="vorontsov.i.e@gmail.com"
ARG BUILD_DATE
ARG BUILD_VERSION

LABEL org.label-schema.schema-version="1.0"
LABEL org.label-schema.build-date=$BUILD_DATE
LABEL org.label-schema.name="PWMBench-PBM"
LABEL org.label-schema.description="PWMBench-PBM benchmark motif performance on PBM data"
LABEL org.label-schema.url="https://github.com/autosome-ru/motif_benchmarks/tree/master/PWMBench-PBM"
LABEL org.label-schema.vcs-url="https://github.com/autosome-ru/motif_benchmarks/tree/master/PWMBench-PBM"
LABEL org.label-schema.version=$BUILD_VERSION

WORKDIR /app
COPY entrypoint.sh  jstacs-2.3.jar  json-simple-1.1.1.jar  /app/
COPY --from=0 /app/*.class /app/
WORKDIR /data
ENTRYPOINT ["/app/entrypoint.sh"]
