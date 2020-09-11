FROM alpine
LABEL maintainer="vorontsov.i.e@gmail.com"
ARG BUILD_DATE
ARG BUILD_VERSION

LABEL org.label-schema.schema-version="1.0"
LABEL org.label-schema.build-date=$BUILD_DATE
LABEL org.label-schema.name="PWMEval-Chip-Peak"
LABEL org.label-schema.description="PWMEval-Chip-Peak benchmark motif performance to recognize ChIP-seq peaks"
LABEL org.label-schema.url="https://github.com/autosome-ru/motif_benchmarks/tree/master/PWMEval-Chip-peak"
LABEL org.label-schema.vcs-url="https://github.com/autosome-ru/motif_benchmarks/tree/master/PWMEval-Chip-peak"
LABEL org.label-schema.version=$BUILD_VERSION

COPY ./chrom_sizes.cpp /source/
RUN apk add --virtual .builddeps --update  alpine-sdk R-dev python bash zlib-dev \
    && mkdir -p /app/ \
     && wget -O /source/pwm_scoring.c  https://raw.githubusercontent.com/gio31415/PWMEval/7973c98a13685fc66984af5dd43048aee7f79456/pwm_scoring.c \
     && gcc -O3 -W -Wall -pedantic -std=gnu99 /source/pwm_scoring.c -o /app/pwm_scoring \
     && g++ -O3 -W -Wall -pedantic /source/chrom_sizes.cpp -o /app/chrom_sizes \
     && rm /source/pwm_scoring.c \
    && mkdir /bedtools && cd /bedtools \
     && wget 'https://github.com/arq5x/bedtools2/releases/download/v2.27.1/bedtools-2.27.1.tar.gz' \
     && tar -zxf bedtools-2.27.1.tar.gz \
     && cd /bedtools/bedtools2 \
     && make \
     && mkdir -p /app && cd /app \
     && cp /bedtools/bedtools2/bin/bedtools /app/bedtools \
     && rm -r /bedtools \
    && apk add R  ttf-ubuntu-font-family rsync \
     && Rscript -e 'install.packages("https://cran.r-project.org/src/contrib/Archive/Rcpp/Rcpp_1.0.1.tar.gz", repos=NULL, INSTALL_opts=c("--no-docs", "--no-help","--no-html"));' \
       && Rscript -e 'install.packages("https://cran.r-project.org/src/contrib/Archive/plyr/plyr_1.8.4.tar.gz", repos=NULL, INSTALL_opts=c("--no-docs", "--no-help","--no-html"));' \
       && Rscript -e 'install.packages("https://cran.r-project.org/src/contrib/Archive/getopt/getopt_1.20.2.tar.gz", repos=NULL, type="source", INSTALL_opts=c("--no-docs", "--no-help","--no-html"));' \
       && Rscript -e 'install.packages("https://cran.r-project.org/src/contrib/Archive/optparse/optparse_1.6.2.tar.gz", repos=NULL, type="source", INSTALL_opts=c("--no-docs", "--no-help","--no-html"));' \
       && Rscript -e 'install.packages("https://cran.r-project.org/src/contrib/bitops_1.0-6.tar.gz", repos=NULL, type="source", INSTALL_opts=c("--no-docs", "--no-help","--no-html"));' \
       && Rscript -e 'install.packages("https://cran.r-project.org/src/contrib/Archive/caTools/caTools_1.17.1.2.tar.gz", repos=NULL, type="source", INSTALL_opts=c("--no-docs", "--no-help","--no-html"));' \
       && Rscript -e 'install.packages("https://cran.r-project.org/src/contrib/rjson_0.2.20.tar.gz", repos=NULL, type="source", INSTALL_opts=c("--no-docs", "--no-help","--no-html"));' \
       && Rscript -e 'install.packages("https://cran.r-project.org/src/contrib/Archive/MASS/MASS_7.3-51.4.tar.gz", repos=NULL, type="source", INSTALL_opts=c("--no-docs", "--no-help","--no-html"));' \
       && Rscript -e 'install.packages("https://cran.r-project.org/src/contrib/Archive/pROC/pROC_1.14.0.tar.gz", repos=NULL, type="source", INSTALL_opts=c("--no-docs", "--no-help","--no-html"));' \
    && apk del .builddeps \
    && rm -rf /var/cache/apk/*
WORKDIR /workdir/
COPY ./calculate_roc.R ./pcm2pfm.R ./utils.R ./motif_preprocessing.R ./peak_preprocessing.R ./assembly_preprocessing.R ./narrowpeak_summit.sh ./peak_center.sh ./download_assembly.sh /app/
ENTRYPOINT ["/app/calculate_roc.R"]
