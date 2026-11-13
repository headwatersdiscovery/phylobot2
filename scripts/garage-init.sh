#!/usr/bin/env bash
# Set up the local Garage node: assign a layout, create the bucket, and import this
# install's access key from .env (scripts/env-secrets.sh generated it). `make up` runs
# it every time, so each step is skipped when already done; a key that changed in .env
# is imported and given the bucket.
set -euo pipefail

set -a; source .env; set +a
garage() { docker compose exec -T garage /garage "$@"; }

if [ -z "${AWS_ACCESS_KEY_ID:-}" ] || [ -z "${AWS_SECRET_ACCESS_KEY:-}" ]; then
  echo "Error: AWS_ACCESS_KEY_ID / AWS_SECRET_ACCESS_KEY are empty in .env; run make env" >&2
  exit 1
fi

# On a first start the container is up before the daemon answers.
for _ in $(seq 60); do
  garage status >/dev/null 2>&1 && break
  sleep 1
done
garage status >/dev/null

if ! garage layout show 2>/dev/null | grep -q "dc1"; then
  node_id=$(garage node id -q | cut -d@ -f1)
  garage layout assign -z dc1 -c 1G "$node_id"
  garage layout apply --version 1
fi

garage bucket info "$S3_BUCKET" >/dev/null 2>&1 || garage bucket create "$S3_BUCKET"

if ! garage key info "$AWS_ACCESS_KEY_ID" >/dev/null 2>&1; then
  garage key import --yes -n phylobot "$AWS_ACCESS_KEY_ID" "$AWS_SECRET_ACCESS_KEY"
  garage bucket allow --read --write --owner "$S3_BUCKET" --key "$AWS_ACCESS_KEY_ID"
  echo "Garage: imported this install's key and gave it bucket '$S3_BUCKET'."
fi
echo "Garage ready: bucket '$S3_BUCKET'."
