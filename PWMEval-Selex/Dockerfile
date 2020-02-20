FROM alpine
COPY ./filter_fasta.cpp /source/
RUN apk add --virtual .builddeps --update  alpine-sdk R-dev \
    && apk add R ttf-ubuntu-font-family \
    && mkdir -p /app/ \
     && wget -O /source/pwm_scoring.c  https://raw.githubusercontent.com/gio31415/PWMEval/8b57b4e60c4cbbd717b38fbb4fae3c9cf770477f/pwm_scoring.c \
     && wget -O /source/seqshuffle.c  https://raw.githubusercontent.com/gio31415/PWMEval/8b57b4e60c4cbbd717b38fbb4fae3c9cf770477f/seqshuffle.c \
     && gcc -O3 -W -Wall -pedantic -std=gnu99 /source/pwm_scoring.c -o /app/pwm_scoring \
     && gcc -O3 -W -Wall -pedantic -std=gnu99 /source/seqshuffle.c -o /app/seqshuffle \
     && g++ -O3 -W -Wall -pedantic /source/filter_fasta.cpp -o /app/filter_fasta \
     && rm /source/pwm_scoring.c /source/seqshuffle.c \
    && Rscript -e 'install.packages(c("rjson", "optparse"), repos = "https://cloud.r-project.org/")' \
    && apk del .builddeps \
    && rm -rf /var/cache/apk/*
WORKDIR /workdir/
COPY ./calculate_roc.R ./fastq2fasta.sh ./filter_redundant.sh ./fasta_multiline_to_singleline.sh ./pcm2pfm.R ./utils.R ./motif_preprocessing.R ./seq_preprocessing.R ./add_flanks.sh  /app/
ENTRYPOINT ["/app/calculate_roc.R"]
