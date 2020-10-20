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
COPY *.rb *.jar download_assembly.sh  /app/
ENTRYPOINT ["ruby", "/app/run_console.rb"]
