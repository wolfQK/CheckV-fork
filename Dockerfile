FROM quay.io/biocontainers/checkv:0.8.1--py_1
LABEL maintainer="Antonio Camargo (antoniop.camargo@lbl.gov), Stephen Nayfach (snayfach@lbl.gov)"
LABEL version="0.8.1"
ENV CHECKVDB="/db/checkv-db-v1.0"
RUN mkdir db
VOLUME [ "/app" ]
WORKDIR /app
RUN checkv download_database /db
ENTRYPOINT [ "checkv" ]