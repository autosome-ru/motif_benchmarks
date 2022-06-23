FROM alpine

COPY filter_fasta.cpp pwm_scoring.c seqshuffle.c  /source/
RUN apk add --virtual .builddeps --update  alpine-sdk R-dev \
    && apk add R ttf-ubuntu-font-family \
    && mkdir -p /app/ \
     && gcc -O3 -W -Wall -pedantic -std=gnu99 /source/pwm_scoring.c -o /app/pwm_scoring \
     && gcc -O3 -W -Wall -pedantic -std=gnu99 /source/seqshuffle.c -o /app/seqshuffle \
     && g++ -O3 -W -Wall -pedantic /source/filter_fasta.cpp -o /app/filter_fasta \
     && rm /source -r \
     && Rscript -e 'install.packages("remotes", repos="http://cran.us.r-project.org");' \
        && Rscript -e 'remotes::install_url("https://cran.r-project.org/src/contrib/Archive/optparse/optparse_1.6.2.tar.gz");' \
        && Rscript -e 'remotes::install_url("https://cran.r-project.org/src/contrib/Archive/rjson/rjson_0.2.0.tar.gz");' \
        && Rscript -e 'remotes::install_url("https://cran.r-project.org/src/contrib/PRROC_1.3.1.tar.gz");' \
    && apk del .builddeps \
    && rm -rf /var/cache/apk/*
RUN cd /app && wget https://github.com/shenwei356/seqkit/releases/download/v0.16.0/seqkit_linux_amd64.tar.gz && tar -zxf seqkit_linux_amd64.tar.gz && rm seqkit_linux_amd64.tar.gz

WORKDIR /workdir/
COPY ./evaluate ./prepare  ./pcm2pfm.R ./pcm2pfm_utils.R ./utils.R ./motif_preprocessing.R ./seq_preprocessing.R ./arglist_options.R ./add_flanks.sh ./roc_pr_curves.R  /app/
ENV PATH="/app:${PATH}"
CMD ["evaluate", "--help"]

ARG BUILD_DATE
ARG BUILD_VERSION

LABEL maintainer="vorontsov.i.e@gmail.com"
LABEL org.label-schema.schema-version="1.0"
LABEL org.label-schema.build-date=$BUILD_DATE
LABEL org.label-schema.name="PWMEval-Selex"
LABEL org.label-schema.description="PWMEval-Selex benchmark motif performance to recognize SELEX-like sequences"
LABEL org.label-schema.url="https://github.com/autosome-ru/motif_benchmarks/tree/master/PWMEval-Selex"
LABEL org.label-schema.vcs-url="https://github.com/autosome-ru/motif_benchmarks/tree/master/PWMEval-Selex"
LABEL org.label-schema.version=$BUILD_VERSION
