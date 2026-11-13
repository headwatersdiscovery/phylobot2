#!/usr/bin/env bash
# Compare .env with .env.example: fail if .env lacks a variable the example defines,
# and mention variables .env has that the example no longer does (safe to delete).
set -euo pipefail

keys() { grep -oE '^[A-Za-z_][A-Za-z0-9_]*=' "$1" | tr -d '=' | sort -u; }

missing=$(comm -23 <(keys .env.example) <(keys .env))
unused=$(comm -13 <(keys .env.example) <(keys .env))

if [ -n "$unused" ]; then
  echo "Note: .env has variables no longer used (safe to delete): $(echo $unused)"
fi
if [ -n "$missing" ]; then
  echo "Error: .env is missing: $(echo $missing)" >&2
  echo "Copy them from .env.example (and set values if you use your own)." >&2
  exit 1
fi
