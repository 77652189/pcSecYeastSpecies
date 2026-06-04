FROM ubuntu:24.04

RUN apt-get update \
    && DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends \
        bash \
        ca-certificates \
        soplex \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /work

