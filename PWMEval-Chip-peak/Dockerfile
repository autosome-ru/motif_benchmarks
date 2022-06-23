FROM alpine

COPY chrom_sizes.cpp pwm_scoring.c  /source/
RUN apk add --virtual .builddeps --update  alpine-sdk R-dev python bash zlib-dev \
    && mkdir -p /app/ \
     && gcc -O3 -W -Wall -pedantic -std=gnu99 /source/pwm_scoring.c -o /app/pwm_scoring \
     && g++ -O3 -W -Wall -pedantic /source/chrom_sizes.cpp -o /app/chrom_sizes \
     && rm -rf /source \
    && mkdir /bedtools && cd /bedtools \
     && wget 'https://github.com/arq5x/bedtools2/releases/download/v2.27.1/bedtools-2.27.1.tar.gz' \
     && tar -zxf bedtools-2.27.1.tar.gz \
     && cd /bedtools/bedtools2 \
     && make \
     && mkdir -p /app && cd /app \
     && cp /bedtools/bedtools2/bin/bedtools /app/bedtools \
     && rm -r /bedtools \
    && apk add R  ttf-ubuntu-font-family rsync \
    && Rscript -e 'install.packages("remotes", repos="http://cran.us.r-project.org");' \
       && Rscript -e 'remotes::install_url("https://cran.r-project.org/src/contrib/Archive/optparse/optparse_1.6.2.tar.gz");' \
       && Rscript -e 'remotes::install_url("https://cran.r-project.org/src/contrib/Archive/rjson/rjson_0.2.0.tar.gz");' \
       && Rscript -e 'remotes::install_url("https://cran.r-project.org/src/contrib/PRROC_1.3.1.tar.gz");' \
    && apk del .builddeps \
    && rm -rf /var/cache/apk/*

WORKDIR /workdir/
COPY ./prepare ./evaluate  ./pcm2pfm.R ./utils.R ./motif_preprocessing.R ./peak_preprocessing.R ./assembly_preprocessing.R ./narrowpeak_summit.sh ./peak_center.sh ./download_assembly.sh ./arglist_options.R ./roc_pr_curves.R  /app/
ENV PATH="/app:${PATH}"
CMD ["evaluate", "--help"]

ARG BUILD_DATE
ARG BUILD_VERSION

LABEL maintainer="vorontsov.i.e@gmail.com"
LABEL org.label-schema.schema-version="1.0"
LABEL org.label-schema.name="PWMEval-Chip-Peak"
LABEL org.label-schema.description="PWMEval-Chip-Peak benchmark motif performance to recognize ChIP-seq peaks"
LABEL org.label-schema.url="https://github.com/autosome-ru/motif_benchmarks/tree/master/PWMEval-Chip-peak"
LABEL org.label-schema.vcs-url="https://github.com/autosome-ru/motif_benchmarks/tree/master/PWMEval-Chip-peak"
LABEL org.label-schema.build-date=$BUILD_DATE
LABEL org.label-schema.version=$BUILD_VERSION
