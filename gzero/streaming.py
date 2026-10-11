"""Large TDMS files in pieces.

A 1 GB TDMS file holds tens of millions of samples per channel. Converting it
in one go needs several GB of memory (every intermediate array is as long as
the file). Here the file is read and converted in chunks instead:

- Loading (:func:`load_streamed`) walks through the file once and keeps only a
  small min/max overview of log G and the piezo for the plot. The recording
  that the program holds is that overview plus the file name and settings.
- Trace detection (:func:`stream_traces`) walks through the file again. Each
  chunk is read with a margin on both sides, converted, and cut into traces
  with the normal code. A chunk keeps only the traces and holds that start in
  its own part of the file, so traces crossing a chunk border are found
  once and complete, as long as the margin is longer than a trace (60 s
  default). The traces and holds are all the analysis tabs need.
- I-V and piezo modulation need the full signals; :func:`read_range` converts
  any part of a file on demand.
"""

from __future__ import annotations

import os
from dataclasses import dataclass, asdict

import numpy as np

from . import bj_io, traces as tr
from .conversion import ConversionSettings, Recording, tdms_to_recording


@dataclass
class StreamSettings:
    threshold_mb: float = 300.0     # TDMS files larger than this are read in pieces
    chunk_s: float = 120.0          # length of one piece (s)
    margin_s: float = 60.0          # read this much before and after each piece (longer than any trace)
    overview_points: int = 400_000  # points of the overview plot

    def to_dict(self) -> dict:
        return asdict(self)


def should_stream(path: str, st: StreamSettings) -> bool:
    return path.lower().endswith(".tdms") and os.path.getsize(path) > st.threshold_mb * 1e6


def _minmax(y: np.ndarray, k: int) -> np.ndarray:
    """Min and max of every block of k samples, in time order (2 points per block)."""
    m = len(y) // k
    if m == 0:
        return y.copy()
    b = y[:m * k].reshape(m, k)
    lo, hi = b.argmin(axis=1), b.argmax(axis=1)
    rows = np.arange(m)
    first = np.where(lo < hi, b[rows, lo], b[rows, hi])
    second = np.where(lo < hi, b[rows, hi], b[rows, lo])
    out = np.column_stack([first, second]).ravel()
    rest = y[m * k:]
    return np.r_[out, [rest.min(), rest.max()]] if len(rest) else out


def chunks(n: int, fs: float, st: StreamSettings):
    """(core start, core stop, read start, read stop) sample ranges covering 0..n."""
    core = max(1, int(st.chunk_s * fs))
    margin = max(0, int(st.margin_s * fs))
    for a in range(0, n, core):
        b = min(n, a + core)
        yield a, b, max(0, a - margin), min(n, b + margin)


def read_range(path: str, cs: ConversionSettings, start: int = 0, stop: int | None = None) -> Recording:
    """Convert part of a TDMS file (all of it with the defaults) to a Recording."""
    with bj_io.TdmsReader(path) as r:
        rec = tdms_to_recording(r.data(start, stop), cs)
    rec.meta["offset"] = int(start)
    return rec


def load_streamed(path: str, cs: ConversionSettings, st: StreamSettings, progress=None) -> Recording:
    """A light Recording of a big TDMS file: the overview of log G (and piezo) plus what is
    needed to read the file again. rec.meta["streamed"] is True; rec.fs is the overview's rate,
    rec.meta["fs_full"] and rec.meta["n_full"] the file's."""
    with bj_io.TdmsReader(path) as r:
        n, fs = r.n, r.fs
        k = max(1, int(np.ceil(n / max(2, st.overview_points // 2))))
        ov_g, ov_p, has_piezo, first_meta, events = [], [], True, None, r.events
        core = max(k, (int(st.chunk_s * fs) // k) * k)          # chunk length: a multiple of k
        for i, a in enumerate(range(0, n, core)):
            b = min(n, a + core)
            rec = tdms_to_recording(r.data(a, b), cs)
            if first_meta is None:
                first_meta = rec.meta
            ov_g.append(_minmax(rec.logG().astype(np.float32), k))
            if rec.piezo is not None:
                ov_p.append(_minmax(np.asarray(rec.piezo, np.float32), k))
            else:
                has_piezo = False
            if progress:
                progress(b, n)
    logG = np.concatenate(ov_g)
    fs_ov = fs * len(logG) / n
    meta = dict(first_meta or {})
    meta.update({"streamed": True, "source": path, "n_full": int(n), "fs_full": float(fs),
                 "conversion": cs.to_dict(), "stream": st.to_dict(),
                 "file_size_mb": os.path.getsize(path) / 1e6})
    return Recording(name=os.path.splitext(os.path.basename(path))[0], fs=fs_ov,
                     G=(10.0 ** logG).astype(np.float32),
                     piezo=np.concatenate(ov_p) if has_piezo and ov_p else None,
                     bias=0.0, current=None, meta=meta, events=events)


def stream_traces(rec: Recording, ss: tr.SegmentSettings, start_index: int = 0, progress=None):
    """Traces and holds of a streamed recording, read chunk by chunk (same result as
    tr.extract_traces on the whole file when the margin is longer than every trace)."""
    path = rec.meta["source"]
    cs = ConversionSettings.from_dict(rec.meta["conversion"])
    st = StreamSettings(**rec.meta["stream"])
    traces, holds, too_long = [], [], 0
    with bj_io.TdmsReader(path) as r:
        n, fs = r.n, r.fs
        for a, b, ra, rb in chunks(n, fs, st):
            part = tdms_to_recording(r.data(ra, rb), cs)
            part.name = rec.name
            t, h = tr.extract_traces(part, ss, start_index=0)
            for x in t + h:
                x.start += ra
                x.stop += ra
            for x in t:
                if a <= x.start < b:
                    if x.stop >= rb - 1 and rb < n:
                        too_long += 1       # cut off by the end of the read window
                    traces.append(x)
            for x in h:
                if a <= x.start < b:
                    holds.append(x)
            if progress:
                progress(b, n)
    for i, x in enumerate(traces):
        x.index = start_index + i
    for i, x in enumerate(holds):
        x.index = i
    rec.meta["traces_cut_by_margin"] = too_long
    return traces, holds
