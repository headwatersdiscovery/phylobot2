"""Nextflow's human-readable durations ('45ms', '12.3s', '1m 3s', '2h 1m') as seconds."""

import re

UNITS = {"ms": 0.001, "s": 1, "m": 60, "h": 3600, "d": 86400}


def seconds(text: str) -> float | None:
    parts = re.findall(r"([\d.]+)\s*(ms|s|m|h|d)\b", text or "")
    if not parts:
        return None
    return sum(float(value) * UNITS[unit] for value, unit in parts)
