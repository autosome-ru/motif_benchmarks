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
RUN mkdir /workdir && chmod ugo+rwx /workdir
WORKDIR /workdir/
COPY *.rb /app/
ENTRYPOINT ["ruby", "/app/run_console.rb"]
