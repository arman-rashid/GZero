"""Writing results into the results folder (txt files + a config json).

I kept the file names of my older analysis files (logHist.txt,
3D_Histogram.txt, ...) so my Origin templates still work.
"""

from __future__ import annotations

import json
import os
from datetime import datetime

import numpy as np


def _jsonable(o):
    if isinstance(o, dict):
        return {str(k): _jsonable(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)):
        return [_jsonable(v) for v in o]
    if isinstance(o, np.ndarray):
        return o.tolist() if o.size < 10000 else f"<array {o.shape}>"
    if isinstance(o, (np.floating, np.integer)):
        return o.item()
    if isinstance(o, float) and not np.isfinite(o):
        return None
    return o


def prepare_folder(base: str, name: str, overwrite: bool = False) -> str:
    """Create base/name; if it exists and not overwrite, append _2, _3, ..."""
    path = os.path.join(base, name)
    if os.path.exists(path) and not overwrite:
        k = 2
        while os.path.exists(f"{path}_{k}"):
            k += 1
        path = f"{path}_{k}"
    os.makedirs(path, exist_ok=True)
    return path


def save_columns(path: str, header: list[str], *cols) -> None:
    arr = np.column_stack([np.asarray(c, float) for c in cols])
    np.savetxt(path, arr, delimiter="\t", header="\t".join(header), comments="", fmt="%.6g")


def save_matrix(path: str, M: np.ndarray) -> None:
    np.savetxt(path, np.asarray(M, float), delimiter="\t", fmt="%.6g")


def save_config(folder: str, prefix: str, config: dict) -> str:
    path = os.path.join(folder, f"{prefix}_config.json")
    config = dict(config)
    config["saved"] = datetime.now().isoformat(timespec="seconds")
    with open(path, "w") as fh:
        json.dump(_jsonable(config), fh, indent=2)
    return path


def save_traces(folder: str, traces, name: str = "Good_Curves") -> str:
    """All selected traces in one tab-separated file: trace index, z (nm), log G."""
    path = os.path.join(folder, f"{name}.txt")
    with open(path, "w") as fh:
        fh.write("trace\tsource\tz_nm\tlogG_G0\n")
        for t in traces:
            for z, y in zip(t.z, t.logG):
                fh.write(f"{t.index}\t{t.rec_name}\t{z:.6g}\t{y:.6g}\n")
    return path


def save_traces_npz(folder: str, traces, name: str = "traces") -> str:
    path = os.path.join(folder, f"{name}.npz")
    np.savez_compressed(
        path,
        index=np.array([t.index for t in traces]),
        source=np.array([t.rec_name for t in traces]),
        kind=np.array([t.kind for t in traces]),
        selected=np.array([t.selected for t in traces]),
        z=np.array([t.z for t in traces], dtype=object),
        logG=np.array([t.logG for t in traces], dtype=object),
    )
    return path


def save_summary(folder: str, lines: list[str], name: str = "Summary.txt") -> str:
    path = os.path.join(folder, name)
    with open(path, "w") as fh:
        fh.write("\n".join(lines) + "\n")
    return path
