#!/usr/bin/env bash
# Remove local pipeline scratch and command-line outputs (make clean).
#
# Removes (all regenerable, none holds submitted jobs' results):
#   work/                  Nextflow's per-task scratch for `make pipeline` / `make golden`
#   results/               outputs of `make pipeline` and golden runs
#   .nextflow/ .nextflow.log*   Nextflow's run history and logs
# Keeps:
#   .nextflow-home/        Nextflow's own framework files (re-downloaded if missing)
#   jobs/                  managed by the worker (kept only for failed web jobs)
#   Docker volumes         Postgres and Garage hold submitted jobs and their results;
#                          never remove them to save space (`docker compose down -v`).
#
# Refuses to run while any pipeline is running: a command-line runner
# (phylobot-runner-run-*), a worker job (phylobot-job-*) or a Nextflow task (nxf-*).
set -euo pipefail

running=$(docker ps --format '{{.Names}}' 2>/dev/null \
  | grep -E '^(phylobot-runner-run-|phylobot-job-|nxf-)' || true)
if [ -n "$running" ]; then
  echo "Not cleaning: a pipeline is running (its scratch is in use):" >&2
  echo "$running" | sed 's/^/  /' >&2
  echo "Wait for it to finish, then run make clean again." >&2
  exit 1
fi

targets=(work results .nextflow)
targets+=(.nextflow.log*)
present=()
for target in "${targets[@]}"; do
  [ -e "$target" ] && present+=("$target")
done
if [ ${#present[@]} -eq 0 ]; then
  echo "Nothing to clean."
  exit 0
fi

freed=$(du -sck "${present[@]}" 2>/dev/null | tail -1 | cut -f1)
rm -rf "${present[@]}"
echo "Removed ${present[*]} ($((freed / 1024)) MB)."
echo "Submitted jobs are unaffected: their results live in Garage, records in Postgres."
