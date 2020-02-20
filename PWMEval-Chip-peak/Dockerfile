FROM alpine
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
    && Rscript -e 'install.packages(c("rjson", "optparse", "MASS", "caTools", "pROC"), repos = "https://cloud.r-project.org/")' \
    && apk del .builddeps \
    && rm -rf /var/cache/apk/*
WORKDIR /workdir/
COPY ./calculate_roc.R ./pcm2pfm.R ./utils.R ./motif_preprocessing.R ./peak_preprocessing.R ./assembly_preprocessing.R ./narrowpeak_summit.sh ./peak_center.sh ./download_assembly.sh /app/
ENTRYPOINT ["/app/calculate_roc.R"]
