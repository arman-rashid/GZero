"""Change-point detection in breaking traces.

Two methods:

Rupture (Chow test), after Hamill et al. 2024 ("Improving single-molecule
conductance measurements with change point detection from the econometrics
toolbox"). Inside the molecular window, log G against distance is fitted by
two straight lines that meet at a break point k instead of one line. The
Chow statistic

    F(k) = [(S_1 - S_2(k)) / p] / [S_2(k) / (n - 2p)],   p = 2

compares the residual sum of squares S_1 of one line with S_2(k) of the two.
The break with the largest F is where the junction ruptures; the plateau
conductance of the trace is the mean log G before it, and the plateau length
the distance from the start of the window to the break (both from the
median, so a few points of the jump into the window don't move them).
A trace counts as a molecular junction if F reaches the threshold and the
line before the break is flat (|slope| below *max plateau slope*, as in the
Junction statistics tab). A plain tunnelling decay can also give a large F,
because its slope changes slightly, but both of its lines are steep, so it
is counted as blank. Because the best of many k is taken, the plain F
distribution would overstate the significance, so the threshold is a value
of F, not a p-value.

Steps (PELT), a multiple-change-point fit in the spirit of the PELT-R idea
(Vaidya, Dyer and Monti): log G is described as a sequence of flat levels.
The Pruned Exact Linear Time algorithm (Killick et al. 2012) finds the
segmentation that minimises

    sum over segments of SSE / sigma^2  +  penalty x 2 ln(n) per change point

with sigma the noise of the trace (from the median absolute difference of
neighbouring points, at least *noise floor*). A tunnelling ramp is cut into
a staircase of short levels by such a fit, so only levels that are flat
themselves (|slope| below *max plateau slope*) and at least *min level
length* long count as plateaus. Every trace then gives its plateau levels
(conductance of each step), their lengths and the step heights between
neighbouring plateaus. Traces are
block-averaged to at most *max points* first, which keeps PELT fast and
doesn't move steps by more than one block.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass, field

import numpy as np

from .traces import Trace


@dataclass
class ChangePointSettings:
    g_high: float = -2.0               # analysis window: log G below this ...
    g_low: float = -5.5                # ... and above this
    min_points: int = 10               # shortest segment (points)
    f_threshold: float = 10.0          # Chow F for a rupture
    max_plateau_slope: float = 2.0     # decades / nm: the part before the break must be flatter
    penalty: float = 1.0               # PELT: multiplies 2 ln(n) per change point
    max_points: int = 400              # PELT: block-average traces to at most this many points
    noise_floor: float = 0.05          # PELT: noise assumed at least this (decades)
    min_step: float = 0.15             # PELT: merge neighbouring levels closer than this (decades)
    min_level_nm: float = 0.05         # PELT: levels shorter than this are not histogrammed
    bins_per_decade: int = 40

    def to_dict(self) -> dict:
        return asdict(self)


@dataclass
class Rupture:
    trace: int
    k: int                 # break index inside the window
    z_break: float         # nm
    F: float
    logG_plateau: float    # median log G before the break
    length: float          # nm from the start of the window to the break
    slope_before: float    # decades / nm
    slope_after: float
    significant: bool


@dataclass
class StepFit:
    trace: int
    bounds: list           # (start, stop) indices in the window, per level
    levels: np.ndarray     # log G of each level
    lengths: np.ndarray    # nm
    slopes: np.ndarray     # decades / nm inside each level
    z: np.ndarray          # distance of the (averaged) points in the window
    y: np.ndarray          # log G of the (averaged) points in the window


@dataclass
class ChangePointResult:
    ruptures: list
    steps: list
    settings: ChangePointSettings
    index: np.ndarray = field(default_factory=lambda: np.zeros(0, int))

    def significant(self):
        return [r for r in self.ruptures if r.significant]

    def summary_lines(self) -> list[str]:
        n = len(self.ruptures)
        if not n:
            return ["No traces with enough points in the window."]
        sig = self.significant()
        out = [f"<b>Rupture (Chow test):</b> {len(sig)} of {n} traces ({100 * len(sig) / n:.1f} %) have a flat "
               f"part (|slope| <= {self.settings.max_plateau_slope:g} dec/nm) that ruptures at F >= "
               f"{self.settings.f_threshold:g}"]
        if sig:
            g = np.array([r.logG_plateau for r in sig])
            out.append(f"plateau log G: median {np.median(g):.3f}, mean {g.mean():.3f} +/- "
                       f"{g.std(ddof=1) / np.sqrt(len(g)) if len(g) > 1 else np.nan:.3f} (SE), spread {g.std():.3f}")
            out.append(f"length to rupture: median {np.median([r.length for r in sig]):.3f} nm")
        if self.steps:
            nl = np.array([int(plateaus(f, self.settings).sum()) for f in self.steps])
            out.append(f"<b>Steps (PELT):</b> {nl.mean():.2f} plateau levels per trace on average; "
                       + ", ".join(f"{k}: {100 * np.mean(nl == k):.0f} %" for k in range(min(4, nl.max() + 1)))
                       + (f", 4 or more: {100 * np.mean(nl >= 4):.0f} %" if nl.max() >= 4 else ""))
        return out


def _cum(x, y):
    return (np.cumsum(np.r_[0.0, np.ones_like(x)]), np.cumsum(np.r_[0.0, x]), np.cumsum(np.r_[0.0, y]),
            np.cumsum(np.r_[0.0, x * x]), np.cumsum(np.r_[0.0, x * y]), np.cumsum(np.r_[0.0, y * y]))


def _sse_line(c, a, b):
    """Residual sum of squares of a straight-line fit to points a..b-1 (arrays of a, b allowed)."""
    n = c[0][b] - c[0][a]
    sx, sy = c[1][b] - c[1][a], c[2][b] - c[2][a]
    sxx, sxy, syy = c[3][b] - c[3][a], c[4][b] - c[4][a], c[5][b] - c[5][a]
    vx = sxx - sx * sx / n
    vy = syy - sy * sy / n
    cxy = sxy - sx * sy / n
    with np.errstate(divide="ignore", invalid="ignore"):
        r = vy - np.where(vx > 0, cxy * cxy / vx, 0.0)
    return np.maximum(r, 0.0)


def window(trace: Trace, s: ChangePointSettings):
    """Distance and log G of the part of a trace inside the window: from the first point
    below g_high to the first point below g_low after it."""
    y = np.asarray(trace.logG, float)
    i0 = np.flatnonzero(y < s.g_high)
    if not len(i0):
        return None, None
    a = int(i0[0])
    below = np.flatnonzero(y[a:] < s.g_low)
    b = a + int(below[0]) if len(below) else len(y)
    return np.asarray(trace.z, float)[a:b], y[a:b]


def chow(z: np.ndarray, y: np.ndarray, min_points: int):
    """Break index with the largest Chow F and that F (two lines against one)."""
    n = len(y)
    m = max(3, int(min_points))
    if n < 2 * m + 1:
        return None, 0.0
    c = _cum(z, y)
    s1 = float(_sse_line(c, 0, n))
    k = np.arange(m, n - m + 1)
    s2 = _sse_line(c, np.zeros_like(k), k) + _sse_line(c, k, np.full_like(k, n))
    p = 2
    with np.errstate(divide="ignore", invalid="ignore"):
        F = ((s1 - s2) / p) / (s2 / (n - 2 * p))
    F = np.nan_to_num(F, nan=0.0, posinf=1e12)
    j = int(np.argmax(F))
    return int(k[j]), float(F[j])


def rupture(trace: Trace, s: ChangePointSettings, k_trace: int = 0) -> Rupture | None:
    z, y = window(trace, s)
    if z is None:
        return None
    k, F = chow(z, y, s.min_points)
    if k is None:
        return None
    sb, sa = _slope(z[:k], y[:k]), _slope(z[k:], y[k:])
    ok = F >= s.f_threshold and abs(sb) <= s.max_plateau_slope
    return Rupture(k_trace, k, float(z[k]), F, float(np.median(y[:k])), float(abs(z[k] - z[0])), sb, sa, bool(ok))


def _slope(z, y) -> float:
    if len(z) < 2 or np.ptp(z) == 0:
        return np.nan
    return float(np.polyfit(z, y, 1)[0])


def block_average(x: np.ndarray, n_max: int):
    m = int(np.ceil(len(x) / max(int(n_max), 1)))
    if m <= 1:
        return np.asarray(x, float)
    n = len(x) // m * m
    return np.asarray(x[:n], float).reshape(-1, m).mean(axis=1)


def noise_sigma(y: np.ndarray) -> float:
    d = np.diff(y)
    if len(d) < 2:
        return 1.0
    s = 1.4826 * np.median(np.abs(d - np.median(d))) / np.sqrt(2)
    return float(s) if s > 0 else float(np.std(y) or 1.0)


def pelt(y: np.ndarray, penalty: float, min_size: int = 2, sigma: float | None = None) -> list[int]:
    """Change points (indices where a new segment starts) of a mean-shift model, by PELT."""
    y = np.asarray(y, float)
    n = len(y)
    if n < 2 * min_size:
        return []
    sg = sigma or noise_sigma(y)
    cs = np.cumsum(np.r_[0.0, y])
    cs2 = np.cumsum(np.r_[0.0, y * y])

    def cost(a, b):          # a: array of starts, b: end
        m = b - a
        return (cs2[b] - cs2[a] - (cs[b] - cs[a]) ** 2 / m) / sg ** 2

    F = np.full(n + 1, np.inf)
    F[0] = -penalty
    last = np.zeros(n + 1, int)
    cand = np.array([0])
    for t in range(1, n + 1):
        young = cand[t - cand < min_size]
        ok = cand[t - cand >= min_size]
        if len(ok):
            c = F[ok] + cost(ok, t)
            v = c + penalty
            j = int(np.argmin(v))
            F[t], last[t] = v[j], ok[j]
            ok = ok[c <= F[t]]                 # PELT pruning
        cand = np.r_[ok, young]
        if np.isfinite(F[t]):
            cand = np.r_[cand, t]
    cps, t = [], n
    while t > 0:
        t = int(last[t])
        if t > 0:
            cps.append(t)
    return sorted(cps)


def steps(trace: Trace, s: ChangePointSettings, k_trace: int = 0) -> StepFit | None:
    z, y = window(trace, s)
    if z is None or len(y) < 2 * s.min_points:
        return None
    zb, yb = block_average(z, s.max_points), block_average(y, s.max_points)
    m = int(np.ceil(len(y) / max(int(s.max_points), 1)))
    min_size = max(2, int(np.ceil(s.min_points / max(m, 1))))
    sg = max(noise_sigma(y) / np.sqrt(max(m, 1)), s.noise_floor)
    cps = pelt(yb, s.penalty * 2 * np.log(len(yb)), min_size, sg)
    edges = [0] + cps + [len(yb)]
    bounds = [(a, b) for a, b in zip(edges[:-1], edges[1:])]
    # merge neighbours closer than min_step
    merged = [bounds[0]]
    for a, b in bounds[1:]:
        pa, pb = merged[-1]
        if abs(yb[a:b].mean() - yb[pa:pb].mean()) < s.min_step:
            merged[-1] = (pa, b)
        else:
            merged.append((a, b))
    levels = np.array([yb[a:b].mean() for a, b in merged])
    dz = np.median(np.abs(np.diff(zb))) if len(zb) > 1 else 0.0
    lengths = np.array([abs(zb[b - 1] - zb[a]) + dz for a, b in merged])
    slopes = np.array([_slope(zb[a:b], yb[a:b]) if b - a > 2 else np.nan for a, b in merged])
    return StepFit(k_trace, merged, levels, lengths, slopes, zb, yb)


def plateaus(f: StepFit, s: ChangePointSettings) -> np.ndarray:
    """Which levels of a step fit are plateaus (flat and long enough)."""
    return (np.abs(np.nan_to_num(f.slopes, nan=np.inf)) <= s.max_plateau_slope) & (f.lengths >= s.min_level_nm)


def analyse(traces: list[Trace], s: ChangePointSettings, do_steps: bool = True) -> ChangePointResult:
    rup, stp = [], []
    for k, t in enumerate(traces):
        r = rupture(t, s, k)
        if r is not None:
            rup.append(r)
        if do_steps:
            f = steps(t, s, k)
            if f is not None:
                stp.append(f)
    return ChangePointResult(rup, stp, s, np.array([t.index for t in traces]))


def level_histogram(res: ChangePointResult):
    """Histogram of the PELT plateau levels weighted by their length (nm), and the step heights
    between neighbouring plateaus of a trace."""
    s = res.settings
    lv = [(l, w) for f in res.steps for l, w, p in zip(f.levels, f.lengths, plateaus(f, s)) if p]
    edges = np.linspace(s.g_low, s.g_high, max(5, int(round((s.g_high - s.g_low) * s.bins_per_decade))) + 1)
    c = 0.5 * (edges[1:] + edges[:-1])
    h = np.histogram([l for l, _ in lv], edges, weights=[w for _, w in lv])[0] if lv else np.zeros(len(c))
    d = []
    for f in res.steps:
        lv_ = f.levels[plateaus(f, s)]
        d += list(np.diff(lv_))
    d = np.array(d)
    return c, h, d
