FROM openjdk:jre-alpine
COPY entrypoint.sh /entrypoint.sh
COPY jstacs-2.3.jar /jstacs-2.3.jar
COPY *.class /
WORKDIR /data
ENTRYPOINT ["/entrypoint.sh"]

# javac --release 8 -cp .:jstacs-2.3.jar PWMBench.java
# docker build -t vorontsovie/pwmeval_pbm -f ~/iogen/projects/greco_benchmark/benchmarks/pbm/Dockerfile ~/iogen/projects/greco_benchmark/benchmarks/pbm/
# docker run --rm --mount type=bind,src=$(pwd),dst=/data vorontsovie/pwmeval_pbm /data/ARX_REF_R1_combinatorial.txt /data/AHR_HUMAN.H11MO.0.B.mat
