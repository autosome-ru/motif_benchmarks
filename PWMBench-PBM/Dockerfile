FROM openjdk:8-jdk-alpine
WORKDIR /app
COPY entrypoint.sh /app/entrypoint.sh
COPY PWMBench.java /app/PWMBench.java
COPY jstacs-2.3.jar /app/jstacs-2.3.jar
RUN javac -cp .:jstacs-2.3.jar PWMBench.java

FROM openjdk:8-jre-alpine
WORKDIR /app
COPY entrypoint.sh /app/entrypoint.sh
COPY jstacs-2.3.jar /app/jstacs-2.3.jar
COPY --from=0 /app/*.class /app/
WORKDIR /data
ENTRYPOINT ["/app/entrypoint.sh"]
