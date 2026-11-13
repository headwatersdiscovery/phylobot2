# syntax=docker/dockerfile:1
FROM python:3.14.7-slim AS base

COPY --from=ghcr.io/astral-sh/uv:0.12.15 /uv /uvx /bin/

ENV UV_COMPILE_BYTECODE=1 \
    UV_LINK_MODE=copy \
    UV_PYTHON_DOWNLOADS=never \
    PYTHONUNBUFFERED=1

WORKDIR /app

# Dependencies first, so code changes don't invalidate this layer.
RUN --mount=type=cache,target=/root/.cache/uv \
    --mount=type=bind,source=uv.lock,target=uv.lock \
    --mount=type=bind,source=pyproject.toml,target=pyproject.toml \
    uv sync --locked --no-install-project --no-dev

COPY . /app
RUN --mount=type=cache,target=/root/.cache/uv \
    uv sync --locked --no-dev

ENV PATH="/app/.venv/bin:$PATH"

# Test runner: the web image plus dev tools, run against the compose Postgres.
FROM base AS test
RUN --mount=type=cache,target=/root/.cache/uv \
    uv sync --locked
CMD ["pytest", "-m", "not golden and not e2e"]

# Job worker: the web image plus the Docker CLI, to launch the Nextflow runner through
# the host Docker socket. Runs as root: socket access is root-equivalent anyway.
FROM base AS worker
COPY --from=docker:29.8.1-cli /usr/local/bin/docker /usr/local/bin/docker
CMD ["python", "manage.py", "run_worker"]

FROM base AS web
RUN useradd --system --uid 1000 --no-create-home phylobot
USER phylobot
EXPOSE 8000
CMD ["gunicorn", "phylobot.wsgi:application", "--bind", "0.0.0.0:8000", "--workers", "2"]
