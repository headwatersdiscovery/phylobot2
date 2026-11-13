# syntax=docker/dockerfile:1
# phylobot-grasp: GRASP (Bodén lab, bnkit) for indel-aware ancestral reconstruction.
#
# bnkit.jar is downloaded at build time and checked against its SHA-256. The bnkit
# repository carries no licence file, so the jar must not be redistributed: build this
# image locally (make grasp-image); never push it to a public registry.
# The wrapper lets Java's heap use 75% of the container's memory (the task's `memory`,
# 2 GB): Java's default of 25% ran out at 80 sequences x 400 residues.
FROM eclipse-temurin:21.0.12.1_1-jre-noble

ARG GRASP_VERSION=2.0.0
ARG GRASP_SHA256=d10453102b6987e2441186b787a9dddad9aebe88dccde37ab2cc1925cca2bc62

RUN mkdir -p /opt/grasp \
    && curl -fsSL -o /opt/grasp/bnkit.jar \
       "https://github.com/bodenlab/bnkit/releases/download/GRASP_v${GRASP_VERSION}/bnkit.jar" \
    && echo "${GRASP_SHA256}  /opt/grasp/bnkit.jar" | sha256sum -c - \
    && printf '#!/bin/sh\nexec java -XX:MaxRAMPercentage=75 -jar /opt/grasp/bnkit.jar "$@"\n' > /usr/local/bin/grasp \
    && chmod +x /usr/local/bin/grasp
