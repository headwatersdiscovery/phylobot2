# syntax=docker/dockerfile:1
# phylobot-tools: the phylobot package (without the web stack) for pipeline tasks.
# Provides the phylobot-validate / phylobot-export / phylobot-summary commands, and
# iqtree3 for phylobot-esr.
# IQ-TREE for ESR (phylobot-esr), whose steps interleave Python and IQ-TREE: the binary
# from the pipeline's own pinned container (pipeline/nextflow.config, pinned to it by a
# test), with the two libraries it loads from ../lib beside it.
FROM quay.io/biocontainers/iqtree:3.1.3--h8471819_0 AS iqtree

FROM python:3.14.7-slim

COPY --from=iqtree /usr/local/bin/iqtree3 /opt/iqtree/bin/iqtree3
COPY --from=iqtree /usr/local/lib/libstdc++.so.6 /usr/local/lib/libgcc_s.so.1 /opt/iqtree/lib/
COPY --from=ghcr.io/astral-sh/uv:0.12.15 /uv /bin/

# Nextflow needs `ps` in every task container to collect task metrics.
RUN apt-get update \
    && apt-get install -y --no-install-recommends procps \
    && rm -rf /var/lib/apt/lists/*

ENV UV_COMPILE_BYTECODE=1 \
    UV_LINK_MODE=copy \
    UV_PYTHON_DOWNLOADS=never \
    UV_PROJECT_ENVIRONMENT=/opt/phylobot

WORKDIR /opt/src

# Runtime dependencies only (no dev or web groups); cached apart from code changes.
RUN --mount=type=cache,target=/root/.cache/uv \
    --mount=type=bind,source=uv.lock,target=uv.lock \
    --mount=type=bind,source=pyproject.toml,target=pyproject.toml \
    uv sync --locked --no-install-project --no-default-groups

COPY pyproject.toml uv.lock README.md LICENSE ./
COPY src ./src
RUN --mount=type=cache,target=/root/.cache/uv \
    uv sync --locked --no-default-groups --no-editable

ENV PATH="/opt/phylobot/bin:/opt/iqtree/bin:$PATH"
WORKDIR /
