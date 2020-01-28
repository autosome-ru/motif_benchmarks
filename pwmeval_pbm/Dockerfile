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

# javac --release 8 -cp .:jstacs-2.3.jar PWMBench.java
# docker build -t vorontsovie/pwmeval_pbm -f ~/iogen/projects/greco_benchmark/benchmarks/pbm/Dockerfile ~/iogen/projects/greco_benchmark/benchmarks/pbm/
# docker run --rm --mount type=bind,src=$(pwd),dst=/data vorontsovie/pwmeval_pbm /data/ARX_REF_R1_combinatorial.txt /data/AHR_HUMAN.H11MO.0.B.mat
