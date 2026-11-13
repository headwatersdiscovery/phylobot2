#!/usr/bin/env bash
# Refuse to (re)start the stack while the worker runs a job (make up, and so make e2e).
#
# `docker compose up --build` recreates the worker whenever its image or config changed,
# which stops every job it is running: the job is marked "interrupted" and has to be
# resubmitted (it cost the owner a job on 2026-10-18). Each running job has a Nextflow
# head container named phylobot-job-<job id>-<attempt> (executors.LocalDockerExecutor); a
# command-line pipeline (make pipeline, make golden) is not affected and not checked.
# FORCE=1 starts the stack anyway.
set -euo pipefail

[ "${FORCE:-}" = "1" ] && exit 0
running=$(docker ps --format '{{.Names}}' 2>/dev/null | grep -E '^phylobot-job-' || true)
if [ -n "$running" ]; then
  echo "Not starting the stack: the worker is running $(echo "$running" | wc -l | tr -d ' ') job(s)," >&2
  echo "and recreating it would stop them (they would have to be resubmitted):" >&2
  echo "$running" | sed 's/^phylobot-job-/  job /' >&2
  echo "Wait for them to finish (My jobs), or run with FORCE=1 to stop them anyway." >&2
  exit 1
fi
