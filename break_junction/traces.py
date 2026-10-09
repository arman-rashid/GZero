"""Cutting a recording into single traces.

With a piezo signal the recording is split where the piezo turns or stops,
without one at the contacts. Whether a piece is breaking or making comes
from the conductance. Each trace is then cut between the high and low
limits and gets a distance axis with zero at the zero-set crossing.

All conductance limits are log10(G/G0).
"""

from __future__ import annotations

from dataclasses import dataclass, asdict

import numpy as np
from scipy.ndimage import uniform_filter1d

from .conversion import Recording

LOG_FLOOR = -10.0


@dataclass
class Trace:
    """One cut trace. z in nm (0 at the zero-set crossing), logG log10(G/G0)."""

    index: int
    kind: str                 # "break" | "make" | "hold"
    rec_name: str
    start: int                # sample range in the recording (after cutting)
    stop: int
    z: np.ndarray
    logG: np.ndarray
    fs: float
    good: bool = True
    reason: str = ""
    manual: bool | None = None   # user override of good/bad (None = automatic)
    cluster: int = -1            # cluster label from the clustering tab (-1 = none)

    @property
    def selected(self) -> bool:
        return self.manual if self.manual is not None else self.good

    @property
    def t(self) -> np.ndarray:
        return np.arange(len(self.logG)) / self.fs

    @property
    def G(self) -> np.ndarray:
        return 10.0 ** self.logG


@dataclass
class SegmentSettings:
    method: str = "piezo"            # "piezo" | "threshold"
    use: str = "break"               # "break" | "make" | "both"
    smooth_ms: float = 10.0          # piezo smoothing for direction detection
    hold_fraction: float = 0.1       # |v| below this fraction of typical speed = hold
    min_segment_ms: float = 30.0     # merge direction runs shorter than this
    classify: str = "conductance"    # "conductance" | "piezo+" | "piezo-" (which direction breaks)
    join_holds: bool = True          # pull-hold-pull -> one trace (holds still kept for noise)
    # distance
    distance_from: str = "piezo"     # "piezo" (smoothed signal) | "index" (points x piezo rate)
    displacement_ratio: float = 1.0  # junction nm per piezo nm (MCBJ attenuation)
    piezo_rate: float = 0.001        # nm per sample for "index" mode / files without piezo
    zero_set: float = -0.3           # log G at which distance = 0 (0.5 G0)
    # cutting limits
    high_check: float = 0.0          # trace must reach above this (contact)
    low_check: float = -4.5          # trace must reach below this
    low_cut: float = -5.0            # cut after first drop below this ...
    additional_points: int = 100     # ... plus this many points
    high_start: float | None = None  # start of the cut trace (None = high_check)
    # quality
    noise_check: bool = False
    noise_points: int = 25           # std of log G over the first points
    noise_limit: float = 0.2         # max std of log G in those points
    max_floor_noise: float | None = None  # optional: max std of log G in last noise_points

    def to_dict(self) -> dict:
        return asdict(self)


# Segmentation

def _runs(state: np.ndarray):
    """Run-length encode an int array -> (starts, stops, values)."""
    if len(state) == 0:
        return np.array([], int), np.array([], int), np.array([], int)
    change = np.flatnonzero(np.diff(state)) + 1
    starts = np.r_[0, change]
    stops = np.r_[change, len(state)]
    return starts, stops, state[starts]


def piezo_segments(piezo: np.ndarray, fs: float, ss: SegmentSettings):
    """Return list of (start, stop, direction) with direction in {+1, -1, 0 (hold)}."""
    w = max(3, int(round(ss.smooth_ms * 1e-3 * fs)) | 1)
    xs = uniform_filter1d(np.asarray(piezo, float), w, mode="nearest")
    xs = uniform_filter1d(xs, w, mode="nearest")
    v = np.gradient(xs) * fs
    moving = np.abs(v)
    scale = np.percentile(moving, 90)
    if scale <= 0:
        return [(0, len(piezo), 0)]
    state = np.where(moving > ss.hold_fraction * scale, np.sign(v), 0).astype(int)

    min_len = max(1, int(ss.min_segment_ms * 1e-3 * fs))
    for _ in range(10):
        starts, stops, vals = _runs(state)
        short = np.flatnonzero((stops - starts) < min_len)
        if len(short) == 0:
            break
        for k in short:
            # absorb short run into the longer neighbour
            left = vals[k - 1] if k > 0 else None
            right = vals[k + 1] if k + 1 < len(vals) else None
            if left is None and right is None:
                continue
            if left is None:
                fill = right
            elif right is None:
                fill = left
            else:
                ll = stops[k - 1] - starts[k - 1]
                rl = stops[k + 1] - starts[k + 1]
                fill = left if ll >= rl else right
            state[starts[k]:stops[k]] = fill
    starts, stops, vals = _runs(state)
    return [(int(a), int(b), int(c)) for a, b, c in zip(starts, stops, vals)]


def threshold_segments(logG: np.ndarray, ss: SegmentSettings, fs: float):
    """Segments for files without piezo: contact runs (logG >= high_check) delimit cycles.

    Between two contact runs the junction breaks then re-forms; split that gap
    at its deepest smoothed point into a break part (+1) and a make part (-1).
    """
    contact = logG >= ss.high_check
    starts, stops, vals = _runs(contact.astype(int))
    contacts = [(a, b) for a, b, v in zip(starts, stops, vals) if v == 1]
    segs = []
    if not contacts:
        return [(0, len(logG), 0)]
    smooth = uniform_filter1d(logG, max(3, int(0.005 * fs)), mode="nearest")
    for (a0, b0), (a1, b1) in zip(contacts[:-1], contacts[1:]):
        gap0, gap1 = b0, a1
        if gap1 - gap0 < 5:
            continue
        mid = gap0 + int(np.argmin(smooth[gap0:gap1]))
        # break segment includes some of the contact before it
        segs.append((max(a0, b0 - (b0 - a0) // 2), mid, +1))
        segs.append((mid, min(b1, a1 + (b1 - a1) // 2), -1))
    # tail after the last contact: a final break
    a_last, b_last = contacts[-1]
    if len(logG) - b_last > 10:
        segs.append((max(a_last, b_last - (b_last - a_last) // 2), len(logG), +1))
    return segs


def _join_holds(segs):
    """Merge move-hold-move with the same direction into one moving segment.

    Hold protocols stop the piezo in the middle of a pull; without joining,
    the breaking trace would end at the hold and miss its tunnelling tail.
    """
    out = []   # (pieces, direction); the hold samples themselves are dropped
    for a, b, d in segs:
        if d != 0 and len(out) >= 2 and out[-1][1] == 0 and out[-2][1] == d:
            out.pop()
            pieces, _ = out.pop()
            out.append((pieces + [(a, b)], d))
        else:
            out.append(([(a, b)], d))
    return out


def _classify(logG_seg: np.ndarray, direction: int, ss: SegmentSettings) -> str:
    if direction == 0:
        return "hold"
    if ss.classify == "piezo+":
        return "break" if direction > 0 else "make"
    if ss.classify == "piezo-":
        return "break" if direction < 0 else "make"
    n = max(3, len(logG_seg) // 20)
    first = np.median(logG_seg[:n])
    last = np.median(logG_seg[-n:])
    if first - last > 1.0:
        return "break"
    if last - first > 1.0:
        return "make"
    return "other"


# Cutting + distance axis

def _cut_breaking(logG: np.ndarray, ss: SegmentSettings):
    """Return (i0, i1, ok, reason) for a breaking trace (conductance decreasing)."""
    if logG.max() < ss.high_check:
        return 0, len(logG), False, "never reaches high limit"
    if logG.min() > ss.low_check:
        return 0, len(logG), False, "never reaches low limit"
    high_start = ss.high_check if ss.high_start is None else ss.high_start
    below = np.flatnonzero(logG < ss.low_cut)
    # first drop below low_cut that happens after the contact
    first_high = int(np.argmax(logG >= high_start))
    below = below[below > first_high]
    i_low = int(below[0]) if len(below) else len(logG) - 1
    above = np.flatnonzero(logG[:i_low + 1] >= high_start)
    if len(above) == 0:
        return 0, len(logG), False, "no contact before break"
    i0 = int(above[-1])
    i1 = min(len(logG), i_low + 1 + ss.additional_points)
    if i1 - i0 < 5:
        return i0, i1, False, "too short"
    return i0, i1, True, ""


def _zero_index(logG: np.ndarray, zero_set: float) -> float:
    """Fractional index of the first downward crossing of zero_set.

    The cut trace starts at the last point above the contact limit, so the
    first crossing is the rupture of the last-atom contact.
    """
    above = logG >= zero_set
    idx = np.flatnonzero(above[:-1] & ~above[1:])
    if len(idx) == 0:
        return 0.0
    i = int(idx[0])
    y0, y1 = logG[i], logG[i + 1]
    frac = (y0 - zero_set) / (y0 - y1) if y0 != y1 else 0.0
    return i + frac


def _distance(piezo_smooth: np.ndarray | None, idx: np.ndarray, ss: SegmentSettings) -> np.ndarray:
    """Distance along the (opening-ordered) sample indices idx."""
    if piezo_smooth is not None and ss.distance_from == "piezo":
        x = piezo_smooth[idx]
        d = x - x[0]
        if len(d) > 1 and d[-1] < 0:
            d = -d
        return d * ss.displacement_ratio
    return np.arange(len(idx)) * ss.piezo_rate


def extract_traces(rec: Recording, ss: SegmentSettings, start_index: int = 0) -> tuple[list[Trace], list[Trace]]:
    """Split rec into traces. Returns (traces, holds).

    traces contains breaking and/or making traces (per ss.use), all
    oriented so that distance increases while the junction opens. holds
    are piezo-hold periods (for conductance-time / noise analysis).
    """
    logG = np.log10(np.clip(rec.G, 10.0 ** LOG_FLOOR, None))
    if ss.method == "piezo" and rec.piezo is not None:
        segs = piezo_segments(rec.piezo, rec.fs, ss)
    else:
        segs = threshold_segments(logG, ss, rec.fs)

    holds = [Trace(i, "hold", rec.name, a, b, np.zeros(b - a), logG[a:b].copy(), rec.fs)
             for i, (a, b, d) in enumerate((s for s in segs if s[2] == 0 and s[1] - s[0] >= 10))]
    groups = _join_holds(segs) if ss.join_holds else [([(a, b)], d) for a, b, d in segs]
    piezo_smooth = None
    if rec.piezo is not None:
        # piezo noise is comparable to the step per sample; smooth for a monotonic distance axis
        w = max(1, int(round(ss.smooth_ms * 1e-3 * rec.fs)) | 1)
        piezo_smooth = uniform_filter1d(np.asarray(rec.piezo, float), w, mode="nearest")

    traces = []
    k = start_index
    for pieces, direction in groups:
        if direction == 0:
            continue
        idx = np.concatenate([np.arange(a, b) for a, b in pieces])
        if len(idx) < 10:
            continue
        kind = _classify(logG[idx], direction, ss)
        if kind == "other":
            continue
        if ss.use != "both" and kind != ss.use:
            continue
        if kind == "make":
            idx = idx[::-1]   # analyse making traces in opening order
        i0, i1, ok, reason = _cut_breaking(logG[idx], ss)
        idx = idx[i0:i1]
        y = logG[idx]
        ra, rb = int(idx.min()), int(idx.max()) + 1
        z = _distance(piezo_smooth, idx, ss)
        zi = _zero_index(y, ss.zero_set)
        z = z - np.interp(zi, np.arange(len(z)), z)
        if ok and ss.noise_check and len(y) > ss.noise_points:
            if np.std(y[:ss.noise_points]) > ss.noise_limit:
                ok, reason = False, "noisy start"
        if ok and ss.max_floor_noise is not None and len(y) > ss.noise_points:
            if np.std(y[-ss.noise_points:]) > ss.max_floor_noise:
                ok, reason = False, "noisy floor"
        traces.append(Trace(k, kind, rec.name, ra, rb, z, y.copy(), rec.fs, ok, reason))
        k += 1
    return traces, holds


def selected(traces: list[Trace]) -> list[Trace]:
    return [t for t in traces if t.selected]


def summary(traces: list[Trace]) -> dict:
    n = len(traces)
    good = sum(t.selected for t in traces)
    reasons = {}
    for t in traces:
        if not t.selected:
            r = t.reason or ("manual" if t.manual is False else "rejected")
            reasons[r] = reasons.get(r, 0) + 1
    return {"total": n, "selected": good, "percentage": 100.0 * good / n if n else 0.0,
            "rejected_by": reasons}
