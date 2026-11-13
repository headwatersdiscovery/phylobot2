#!/usr/bin/env bash
# Give this install its own credentials: fill each secret that is empty in .env with a
# random value. .env.example ships them empty, so no two installs share a credential and
# nothing in the repository works as a key anywhere. A value already set is never
# changed (Postgres and Garage keep the credentials they were first set up with).
set -euo pipefail

ENV_FILE=${1:-.env}

hex() { od -An -tx1 -N "$1" /dev/urandom | tr -d ' \n'; }

# name and generator: Garage wants key IDs "GK" + 24 hex and 64-hex secrets.
generate() {
  case "$1" in
    DJANGO_SECRET_KEY) hex 32 ;;
    POSTGRES_PASSWORD) hex 16 ;;
    AWS_ACCESS_KEY_ID) printf 'GK%s' "$(hex 12)" ;;
    AWS_SECRET_ACCESS_KEY | GARAGE_RPC_SECRET | GARAGE_ADMIN_TOKEN) hex 32 ;;
  esac
}

names=(DJANGO_SECRET_KEY POSTGRES_PASSWORD)
# Storage keys only for the local Garage: on AWS, AWS_* are left empty on purpose so
# boto3 uses the machine's IAM role, and a generated key would break that.
if grep -qE '^S3_ENDPOINT_URL=http://garage:' "$ENV_FILE"; then
  names+=(AWS_ACCESS_KEY_ID AWS_SECRET_ACCESS_KEY GARAGE_RPC_SECRET GARAGE_ADMIN_TOKEN)
fi

filled=()
for name in "${names[@]}"; do
  # Only a line that is present and empty ("NAME=") is filled; a missing line is
  # check-env.sh's to report.
  grep -qE "^${name}=$" "$ENV_FILE" || continue
  value=$(generate "$name")
  # awk, not sed -i: GNU and BSD sed disagree on -i.
  tmp=$(mktemp "${ENV_FILE}.XXXXXX")
  awk -v n="$name" -v v="$value" 'BEGIN { FS = OFS = "=" } $0 == n "=" { print n, v; next } { print }' \
    "$ENV_FILE" >"$tmp"
  chmod 600 "$tmp"
  mv "$tmp" "$ENV_FILE"
  filled+=("$name")
done

if [ ${#filled[@]} -gt 0 ]; then
  echo "Generated credentials for this install in $ENV_FILE: ${filled[*]}"
fi
