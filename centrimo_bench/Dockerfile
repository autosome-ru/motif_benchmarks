FROM vorontsovie/meme:1.0
RUN apt-get update && \
    apt-get install --yes ruby2.7 wget && \
    apt-get install --yes build-essential && \
    mkdir /bedtools && cd /bedtools && \
    wget 'https://github.com/arq5x/bedtools2/releases/download/v2.27.1/bedtools-2.27.1.tar.gz' && \
    tar -zxf bedtools-2.27.1.tar.gz && \
    cd /bedtools/bedtools2 && \
    make && \
    mkdir -p /app && cd /app && \
    cp /bedtools/bedtools2/bin/bedtools /app/bedtools && \
    rm -r /bedtools && \
    apt-get remove --yes build-essential && \
    apt-get autoremove --yes && \
    rm -rf /var/lib/apt/lists/*

COPY *.rb prepare evaluate /app/
ENV PATH="/app:${PATH}"
CMD ["evaluate", "--help"]

ARG BUILD_DATE
ARG BUILD_VERSION

LABEL maintainer="vorontsov.i.e@gmail.com"
LABEL org.label-schema.schema-version="1.0"
LABEL org.label-schema.name="centrimo_bench"
LABEL org.label-schema.description="centrimo_bench benchmark motif performance to recognize ChIP-seq peaks"
LABEL org.label-schema.url="https://github.com/autosome-ru/motif_benchmarks/tree/master/centrimo_bench"
LABEL org.label-schema.vcs-url="https://github.com/autosome-ru/motif_benchmarks/tree/master/centrimo_bench"
LABEL org.label-schema.build-date=$BUILD_DATE
LABEL org.label-schema.version=$BUILD_VERSION
