"""Recorded Open Tree of Life answers (record.py) and the transport that replays them."""

import json
from pathlib import Path
from typing import Any

from phylobot.species.opentree import Answer

FIXTURES = Path(__file__).parent


class Replay:
    """A transport answering from a recorded scenario; any other request fails."""

    def __init__(self, name: str):
        self.exchanges = json.loads((FIXTURES / f"{name}.json").read_text())
        self.calls: list[str] = []

    def __call__(self, endpoint: str, payload: dict[str, Any]) -> Answer:
        self.calls.append(endpoint)
        for exchange in self.exchanges:
            if exchange["endpoint"] == endpoint and exchange["payload"] == payload:
                return Answer(exchange["status"], exchange["body"])
        raise AssertionError(f"unrecorded request {endpoint} {payload}")
