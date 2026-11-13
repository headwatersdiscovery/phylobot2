"""Rigid superposition of one structure on another from paired C-alpha atoms.

The pairs come from the job's alignment: two residues, one in each structure, that
belong to the same alignment column. So two different proteins of the family, or an
extant structure and a model of an ancestor, are laid on each other by homology, as
PyMOL's `align` does from a sequence alignment: fit (Kabsch), drop the pairs that are far
apart (a domain that has moved, a loop), fit again, until the pairs kept no longer change.

The page uses it to put its panels in one frame, so that one camera serves them all.
"""

from dataclasses import dataclass

import numpy as np

MIN_PAIRS = 4
ROUNDS = 8
REJECT = 2.0  # pairs further apart than this many RMSDs are left out of the next fit
KEEP_WITHIN = 3.0  # ... but never pairs closer than this (Angstrom)


@dataclass(frozen=True)
class Fit:
    rotation: tuple[tuple[float, float, float], ...]  # 3 x 3, applied as R @ x
    translation: tuple[float, float, float]
    rmsd: float  # over the pairs kept
    kept: int
    pairs: int

    @property
    def matrix(self) -> list[float]:
        """Twelve numbers for structure.js: the rotation by rows, then the translation."""
        flat = [value for row in self.rotation for value in row] + list(self.translation)
        return [round(value, 5) for value in flat]


def kabsch(moving: np.ndarray, fixed: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """The rotation R and translation t for which R @ moving + t is closest to `fixed`."""
    centre_m, centre_f = moving.mean(axis=0), fixed.mean(axis=0)
    covariance = (moving - centre_m).T @ (fixed - centre_f)
    u, _, vt = np.linalg.svd(covariance)
    sign = np.sign(np.linalg.det(vt.T @ u.T)) or 1.0  # never a mirror image
    rotation = vt.T @ np.diag([1.0, 1.0, sign]) @ u.T
    return rotation, centre_f - rotation @ centre_m


def fit(
    moving: list[tuple[float, float, float]], fixed: list[tuple[float, float, float]]
) -> Fit | None:
    """Lay `moving` on `fixed` (paired points); None with too few pairs to mean anything."""
    if len(moving) != len(fixed) or len(moving) < MIN_PAIRS:
        return None
    a, b = np.array(moving, dtype=float), np.array(fixed, dtype=float)
    keep = np.ones(len(a), dtype=bool)
    rotation, translation, rmsd = np.eye(3), np.zeros(3), 0.0
    for _ in range(ROUNDS):
        rotation, translation = kabsch(a[keep], b[keep])
        distance = np.linalg.norm(a @ rotation.T + translation - b, axis=1)
        rmsd = float(np.sqrt((distance[keep] ** 2).mean()))
        # Far is judged two ways and the stricter one counts: against the RMSD (PyMOL's
        # rule), and against the middle of the distances, which a displaced fifth of the
        # pairs cannot drag along with it as it does the RMSD.
        middle = float(np.median(distance[keep]))
        spread = float(np.median(np.abs(distance[keep] - middle))) * 1.4826
        closer = distance <= max(min(REJECT * rmsd, middle + 3.0 * spread), KEEP_WITHIN)
        if closer.sum() < MIN_PAIRS or (closer == keep).all():
            break
        keep = closer
    if not np.isfinite(rotation).all() or not np.isfinite(translation).all():
        return None
    return Fit(
        tuple((float(row[0]), float(row[1]), float(row[2])) for row in rotation),
        (float(translation[0]), float(translation[1]), float(translation[2])),
        rmsd,
        int(keep.sum()),
        len(a),
    )
