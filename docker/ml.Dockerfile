# syntax=docker/dockerfile:1
# phylobot-ml: phylobot-tools plus the neural property predictors:
# metapredict (disorder, MIT) and ESM-2 (fair-esm, MIT), on CPU-only torch. A superset
# of phylobot-tools, so the PROPERTIES step runs here whenever it is asked for one of
# these properties, whatever else it measures in the same task.
#
# The ESM-2 weights are fetched at build time and checked against their SHA-256, so a
# pipeline task never downloads anything (fair-esm would otherwise fetch them through
# torch.hub on first use). Meta publishes them under MIT, like the code.
FROM python:3.14.7-slim AS build

COPY --from=ghcr.io/astral-sh/uv:0.12.15 /uv /bin/

# metapredict publishes Linux sources only (its wheels are macOS), so it compiles here.
RUN apt-get update \
    && apt-get install -y --no-install-recommends build-essential \
    && rm -rf /var/lib/apt/lists/*

ENV UV_COMPILE_BYTECODE=1 \
    UV_LINK_MODE=copy \
    UV_PYTHON_DOWNLOADS=never \
    UV_PROJECT_ENVIRONMENT=/opt/phylobot

WORKDIR /opt/src

RUN --mount=type=cache,target=/root/.cache/uv \
    --mount=type=bind,source=uv.lock,target=uv.lock \
    --mount=type=bind,source=pyproject.toml,target=pyproject.toml \
    uv sync --locked --no-install-project --no-default-groups --group ml

COPY pyproject.toml uv.lock README.md LICENSE ./
COPY src ./src
RUN --mount=type=cache,target=/root/.cache/uv \
    uv sync --locked --no-default-groups --group ml --no-editable

# ESM-2 35M (Lin et al. 2023), the model phylobot.properties_ml uses; the regression
# file is fair-esm's contact head, which its loader insists on finding beside the model.
RUN mkdir -p /opt/esm && python - <<'EOF'
import hashlib, urllib.request
BASE = "https://dl.fbaipublicfiles.com/fair-esm"
FILES = {
    "models/esm2_t12_35M_UR50D.pt":
        "7f21e80e61d16a71735163ef555d3009afb0c98da74c48e29df08606973cc55e",
    "regression/esm2_t12_35M_UR50D-contact-regression.pt":
        "16641e05d830d0ce863dd152dbb8c2f3ddfa3c3ec2a66080152c8abad01d8585",
}
for path, expected in FILES.items():
    data = urllib.request.urlopen(f"{BASE}/{path}").read()
    got = hashlib.sha256(data).hexdigest()
    if got != expected:
        raise SystemExit(f"{path}: SHA-256 {got}, expected {expected}")
    with open("/opt/esm/" + path.split("/")[-1], "wb") as handle:
        handle.write(data)
EOF

# The predictors' own tests (tests/test_properties_ml.py), which need torch and the
# weights and so cannot run in the web image's test service: `make ml-test`. Built only
# when asked for (--target ml-test), never part of the image the pipeline runs.
FROM build AS ml-test
RUN --mount=type=cache,target=/root/.cache/uv \
    uv sync --locked --no-default-groups --group ml --group dev --no-editable
ENV PATH="/opt/phylobot/bin:$PATH" \
    PHYLOBOT_ESM_MODEL=/opt/esm/esm2_t12_35M_UR50D.pt \
    OMP_NUM_THREADS=1 \
    MKL_NUM_THREADS=1 \
    MPLCONFIGDIR=/tmp/matplotlib

# The image the pipeline runs (the last stage, so a plain build produces it).
FROM python:3.14.7-slim

# Nextflow needs `ps` in every task container to collect task metrics.
RUN apt-get update \
    && apt-get install -y --no-install-recommends procps \
    && rm -rf /var/lib/apt/lists/*

COPY --from=build /opt/phylobot /opt/phylobot
COPY --from=build /opt/esm /opt/esm

# One thread: tasks are single-CPU (nextflow.config) and, measured under amd64 emulation
# on Apple Silicon, more torch threads made ESM-2 no faster and metapredict 40x slower.
ENV PATH="/opt/phylobot/bin:$PATH" \
    PHYLOBOT_ESM_MODEL=/opt/esm/esm2_t12_35M_UR50D.pt \
    OMP_NUM_THREADS=1 \
    MKL_NUM_THREADS=1 \
    MPLCONFIGDIR=/tmp/matplotlib
WORKDIR /
