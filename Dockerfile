FROM quay.io/biocontainers/checkv:0.7.0--py_1
LABEL maintainer="Stephen Nayfach, snayfach@lbl.gov"
LABEL version="1.0"
ENV CHECKVDB="/db/checkv-db-v0.6"
RUN mkdir db
VOLUME [ "/app" ]
WORKDIR /app
RUN checkv download_database /db
ENTRYPOINT [ "checkv" ]