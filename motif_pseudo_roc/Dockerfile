FROM ruby:2.6-alpine
COPY ./chrom_sizes.cpp /source/
RUN apk add --update openjdk8-jre rsync \
	&& apk add --virtual .builddeps --update  alpine-sdk ruby-dev bash python2 \
	&& mkdir -p /app/ \
	 && g++ -O3 -W -Wall -pedantic /source/chrom_sizes.cpp -o /app/chrom_sizes \
	&& gem install json --no-document \
	&& mkdir /bedtools && cd /bedtools \
     && wget 'https://github.com/arq5x/bedtools2/releases/download/v2.27.1/bedtools-2.27.1.tar.gz' \
     && tar -zxf bedtools-2.27.1.tar.gz \
     && cd /bedtools/bedtools2 \
     && make \
     && mkdir -p /app && cd /app \
     && cp /bedtools/bedtools2/bin/bedtools /app/bedtools \
     && rm -r /bedtools \
	&& apk del .builddeps

WORKDIR /workdir/
COPY *.rb *.jar download_assembly.sh prepare evaluate /app/
ENV PATH="/app:${PATH}"
CMD ["evaluate", "--help"]

ARG BUILD_DATE
ARG BUILD_VERSION

LABEL maintainer="vorontsov.i.e@gmail.com"
LABEL org.label-schema.schema-version="1.0"
LABEL org.label-schema.name="motif_pseudo_roc"
LABEL org.label-schema.description="motif_pseudo_roc benchmark motif performance to recognize ChIP-seq peaks"
LABEL org.label-schema.url="https://github.com/autosome-ru/motif_benchmarks/tree/master/motif_pseudo_roc"
LABEL org.label-schema.vcs-url="https://github.com/autosome-ru/motif_benchmarks/tree/master/motif_pseudo_roc"
LABEL org.label-schema.build-date=$BUILD_DATE
LABEL org.label-schema.version=$BUILD_VERSION
