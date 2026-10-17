"""Blinking: molecules jumping in and out of a fixed gap (I(t) method).

With the tip held at a fixed height above the surface (STM) or a fixed
electrode gap (MCBJ), the tunnelling current is constant until a molecule
bridges the gap. The current then jumps up to the molecular level and drops
back when the bond breaks. This "blinking" (also called the I(t) or
telegraph method, Haiss et al. 2003, 2006) gives the single-molecule
conductance without pulling, and the on and off times give the lifetime
of the junction and the rate at which it forms.

The analysis, for every segment (a hold, or a whole recording):

1. log G is smoothed with a moving average over *smooth* ms.
2. The baseline (the tunnelling current without a molecule) is a low
   percentile of log G in blocks of *baseline window* s, interpolated in
   between, so it follows slow drift of the gap and is not pulled up by
   the blinks themselves.
3. A blink starts when log G rises more than *up threshold* decades above
   the baseline and ends when it falls below *down threshold* (hysteresis,
   so noise on a plateau doesn't split it). Blinks shorter than *min on*
   are dropped and gaps shorter than *min off* are closed.
4. For every blink: start, duration, its conductance (median log G while
   on), the baseline under it and the jump height. Blinks whose level is
   outside the molecular window are dropped (metal contact, noise).

From all blinks: the blink conductance histogram with a Gaussian fit (the
molecular conductance), the on-time and off-time distributions with
exponential fits (lifetimes tau_on and tau_off, maximum likelihood with the
shortest accepted time as cut-off), the blink rate and the fraction of time
a molecule is in the gap, and a histogram of log G aligned on the blink
starts.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass, field

import numpy as np
from scipy.optimize import curve_fit


@dataclass
class BlinkSettings:
    source: str = "recording"          # "recording" (whole file / from-to) | "holds"
    smooth_ms: float = 1.0             # moving average before the detection
    baseline_s: float = 0.5            # block length for the baseline
    baseline_percentile: float = 20.0  # low percentile of each block = baseline
    up_threshold: float = 0.6          # decades above the baseline: blink starts
    down_threshold: float = 0.3        # decades above the baseline: blink ends
    min_on_ms: float = 1.0
    min_off_ms: float = 0.5
    g_min: float = -6.0                # molecular window for the blink level (log G/G0)
    g_max: float = -1.0
    bins_per_decade: int = 40          # blink conductance histogram
    time_bins: int = 40                # on / off time histograms
    align_before_ms: float = 10.0      # blink-aligned histogram
    align_after_ms: float = 50.0

    def to_dict(self) -> dict:
        return asdict(self)


@dataclass
class Blink:
    segment: int
    start: int            # sample
    stop: int
    duration_s: float
    logG: float           # median log G while on
    baseline: float       # baseline log G at the blink
    height: float         # logG - baseline (decades)


@dataclass
class BlinkResult:
    blinks: list
    off_times: np.ndarray         # s, between blinks inside one segment (censored ones left out)
    total_time: float             # s analysed
    n_segments: int
    rejected: int                 # blinks outside the molecular window
    settings: BlinkSettings
    peak: tuple = (np.nan, np.nan, np.nan)      # Gaussian centre, its error, sigma
    tau_on: tuple = (np.nan, np.nan)
    tau_off: tuple = (np.nan, np.nan)
    hist: tuple = ()                            # (centres, counts) of the blink levels
    aligned: tuple = ()                         # (t centres, g centres, H)
    segment_names: list = field(default_factory=list)

    @property
    def on_times(self) -> np.ndarray:
        return np.array([b.duration_s for b in self.blinks])

    @property
    def levels(self) -> np.ndarray:
        return np.array([b.logG for b in self.blinks])

    def summary_lines(self) -> list[str]:
        n = len(self.blinks)
        if not n:
            return [f"No blinks in {self.total_time:.1f} s ({self.n_segments} segments). Lower the up threshold or "
                    "check the molecular window."]
        on = self.on_times.sum()
        m, me, sg = self.peak
        out = [f"<b>{n} blinks</b> in {self.total_time:.2f} s ({self.n_segments} segments): "
               f"{n / max(self.total_time, 1e-12):.2f} per s, molecule in the gap {100 * on / max(self.total_time, 1e-12):.1f} % "
               f"of the time; {self.rejected} outside log G [{self.settings.g_min:g}, {self.settings.g_max:g}] dropped"]
        if np.isfinite(m):
            out.append(f"<b>blink conductance log G = {m:.3f} +/- {me:.3f}</b> (sigma {sg:.3f}), "
                       f"G = {10 ** m:.3g} G0")
        out.append(f"median jump height {np.median([b.height for b in self.blinks]):.2f} decades")
        if np.isfinite(self.tau_on[0]):
            out.append(f"on time (lifetime): tau_on = {1e3 * self.tau_on[0]:.3g} +/- {1e3 * self.tau_on[1]:.2g} ms, "
                       f"median {1e3 * np.median(self.on_times):.3g} ms")
        if np.isfinite(self.tau_off[0]):
            out.append(f"off time: tau_off = {1e3 * self.tau_off[0]:.3g} +/- {1e3 * self.tau_off[1]:.2g} ms "
                       f"({len(self.off_times)} gaps)")
        return out


def moving_average(y: np.ndarray, n: int) -> np.ndarray:
    n = int(n)
    if n <= 1 or len(y) < n:
        return np.asarray(y, float)
    c = np.cumsum(np.r_[0.0, np.asarray(y, float)])
    out = (c[n:] - c[:-n]) / n
    pad = n - 1
    return np.r_[np.full(pad // 2, out[0]), out, np.full(pad - pad // 2, out[-1])]


def baseline(logG: np.ndarray, fs: float, s: BlinkSettings) -> np.ndarray:
    """Low percentile of log G in blocks of baseline_s, linearly interpolated."""
    n = len(logG)
    b = int(max(10, s.baseline_s * fs))
    if b >= n:
        return np.full(n, np.percentile(logG, s.baseline_percentile))
    starts = np.arange(0, n, b)
    vals = np.array([np.percentile(logG[a:a + b], s.baseline_percentile) for a in starts])
    centres = np.minimum(starts + b / 2, n - 1)
    return np.interp(np.arange(n), centres, vals)


def _runs(mask: np.ndarray):
    """(start, stop) of the runs of True in mask."""
    d = np.diff(np.r_[0, mask.astype(np.int8), 0])
    return np.flatnonzero(d == 1), np.flatnonzero(d == -1)


def detect(logG: np.ndarray, fs: float, s: BlinkSettings, segment: int = 0):
    """Blinks in one segment: (blinks, smoothed log G, baseline, rejected count)."""
    y = moving_average(np.asarray(logG, float), max(1, int(round(s.smooth_ms * 1e-3 * fs))))
    base = baseline(y, fs, s)
    h = y - base
    # hysteresis: on from a crossing of up until h drops below down
    above_down = h > s.down_threshold
    st, sp = _runs(above_down)
    keep = [(a, b) for a, b in zip(st, sp) if np.any(h[a:b] > s.up_threshold)]
    # close short gaps
    min_off = max(1, int(round(s.min_off_ms * 1e-3 * fs)))
    merged = []
    for a, b in keep:
        if merged and a - merged[-1][1] < min_off:
            merged[-1] = (merged[-1][0], b)
        else:
            merged.append((a, b))
    min_on = max(1, int(round(s.min_on_ms * 1e-3 * fs)))
    out, rejected = [], 0
    for a, b in merged:
        if b - a < min_on:
            continue
        if a == 0 or b == len(y):         # cut by the start / end of the segment: length unknown
            continue
        lev = float(np.median(y[a:b]))
        if not s.g_min <= lev <= s.g_max:
            rejected += 1
            continue
        bl = float(np.median(base[a:b]))
        out.append(Blink(segment, int(a), int(b), (b - a) / fs, lev, bl, lev - bl))
    return out, y, base, rejected


def exp_lifetime(t: np.ndarray):
    """Maximum-likelihood lifetime of an exponential above the shortest time: (tau, error)."""
    t = np.asarray(t, float)
    if len(t) < 3:
        return np.nan, np.nan
    tau = float(np.mean(t - t.min())) * len(t) / max(len(t) - 1, 1)
    return tau, tau / np.sqrt(len(t))


def _gauss(x, a, m, sg):
    return a * np.exp(-(x - m) ** 2 / (2 * sg ** 2))


def gauss_peak(x: np.ndarray, lo: float, hi: float, per_decade: int):
    """Histogram of x and a Gaussian fit of its main peak: ((centre, error, sigma), (centres, counts))."""
    edges = np.linspace(lo, hi, max(5, int(round((hi - lo) * per_decade))) + 1)
    h, _ = np.histogram(x, edges)
    c = 0.5 * (edges[1:] + edges[:-1])
    if len(x) < 3 or h.max() == 0:
        return (np.nan, np.nan, np.nan), (c, h)
    i = int(np.argmax(h))
    sg0 = max(np.std(x), edges[1] - edges[0])
    try:
        popt, pcov = curve_fit(_gauss, c, h, p0=[h[i], c[i], sg0], maxfev=4000)
        m, sg = float(popt[1]), abs(float(popt[2]))
        if not lo <= m <= hi:
            raise RuntimeError
        sel = np.abs(x - m) < 2 * sg
        err = sg / np.sqrt(max(sel.sum(), 1))
        return (m, err, sg), (c, h)
    except (RuntimeError, ValueError):
        return (float(np.median(x)), float(np.std(x) / np.sqrt(len(x))), float(np.std(x))), (c, h)


def aligned_histogram(segments, blinks, s: BlinkSettings, bins_per_decade: int = 30):
    """2D histogram of log G against time from the blink start (counts per blink)."""
    if not blinks:
        return ()
    fs = segments[0][1]
    nb = int(round(s.align_before_ms * 1e-3 * fs))
    na = int(round(s.align_after_ms * 1e-3 * fs))
    t = (np.arange(-nb, na) / fs) * 1e3
    lo, hi = s.g_min - 1.0, s.g_max + 0.5
    ge = np.linspace(lo, hi, int(round((hi - lo) * bins_per_decade)) + 1)
    te = np.r_[t - 0.5e3 / fs, t[-1] + 0.5e3 / fs]
    H = np.zeros((len(t), len(ge) - 1))
    for b in blinks:
        y = segments[b.segment][0]
        a, e = b.start - nb, b.start + na
        idx = np.arange(a, e)
        ok = (idx >= 0) & (idx < len(y))
        H += np.histogram2d(t[ok], y[idx[ok]], [te, ge])[0]
    H /= len(blinks)
    return t, 0.5 * (ge[1:] + ge[:-1]), H


def analyse(segments, s: BlinkSettings, names=None) -> BlinkResult:
    """segments: list of (logG array, fs)."""
    blinks, offs, rej, total = [], [], 0, 0.0
    smoothed = []
    for k, (lg, fs) in enumerate(segments):
        b, y, _, r = detect(lg, fs, s, k)
        smoothed.append((y, fs))
        rej += r
        total += len(lg) / fs
        blinks += b
        offs += [(b2.start - b1.stop) / fs for b1, b2 in zip(b[:-1], b[1:])]
    res = BlinkResult(blinks, np.array(offs), total, len(segments), rej, s, segment_names=list(names or []))
    if blinks:
        res.peak, res.hist = gauss_peak(res.levels, s.g_min, s.g_max, s.bins_per_decade)
        res.tau_on = exp_lifetime(res.on_times)
        res.tau_off = exp_lifetime(res.off_times)
        res.aligned = aligned_histogram(smoothed, blinks, s)
    return res


def simulate(duration_s: float = 20.0, fs: float = 10000.0, logG_mol: float = -3.5, sigma_mol: float = 0.15,
             logG_base: float = -5.2, tau_on: float = 0.02, tau_off: float = 0.1, noise: float = 0.08,
             drift: float = 0.3, seed: int = 0):
    """Telegraph-like log G(t) with exponential on / off times and a slowly drifting baseline.
    Returns (logG, fs, list of true (start, stop, level))."""
    rng = np.random.default_rng(seed)
    n = int(duration_s * fs)
    t = np.arange(n) / fs
    base = logG_base + drift * np.sin(2 * np.pi * t / duration_s)
    y = base.copy()
    truth = []
    i = int(rng.exponential(tau_off) * fs)
    while i < n:
        d = max(1, int(rng.exponential(tau_on) * fs))
        lev = rng.normal(logG_mol, sigma_mol)
        if i + d < n:
            y[i:i + d] = lev
            truth.append((i, i + d, lev))
        i += d + max(1, int(rng.exponential(tau_off) * fs))
    y += rng.normal(0, noise, n)
    return y, fs, truth
