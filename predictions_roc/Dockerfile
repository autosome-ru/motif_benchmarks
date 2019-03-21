FROM ruby:2.6-alpine
RUN apk add --virtual .builddeps --update  alpine-sdk ruby-dev && \
	gem install json --no-document && \
	apk del .builddeps
WORKDIR /workdir/
COPY .  /app/
ENTRYPOINT ["ruby", "/app/run.rb"]
# CMD ["/workdir/config.json", "/workdir/persistent/result.json"]
