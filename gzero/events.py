"""Flickering and mechanical events in conductance-time traces.

Follows Rashid et al., JACS 147, 830 (2025), SI section 2.6:

1. G(t) (in G0) goes through a Butterworth low-pass (200 Hz, order 2) and
   then a Savitzky-Golay filter (polynomial order 20, 20 side points).
2. The filtered trace is differentiated with the second-order central
   difference, dG/dt = (G[i+1] - G[i-1]) / (2 dt). dt is a setting (20 ms in
   the paper, which is how the thresholds below are scaled).
3. Peaks of the derivative are found as LabVIEW's Peak Detector does it: a
   quadratic is fitted over "width" points (50) around every sample and a
   peak counts if the fitted maximum reaches the threshold. Above the
   mechanical threshold (0.002) it is a mechanical event, between the
   flicker threshold (0.0007) and that a flicker. By default only rises
   (positive peaks) are counted, as the Peak Detector does; "both" also
   counts drops. "minimum width" instead keeps peaks whose width at half
   height is at least "width" samples (scipy's find_peaks).
4. Events where the junction is outside the molecular conductance range
   (metal contact, noise floor) are false positives and dropped.

A trace is "mechanical" if it has a mechanical event, "flickering" if it has
flickers only, and "quiet" otherwise.

The module also makes the 2D conductance-time histogram of the holds with
the most probable conductance at each time (Gaussian fit per time bin),
used in the paper for the mechanical modulation (SI Fig. S23).
"""

from __future__ import annotations

from dataclasses import dataclass, asdict, field

import numpy as np
from scipy.optimize import curve_fit
from scipy.signal import butter, find_peaks, sosfiltfilt

from .traces import Trace


@dataclass
class EventSettings:
    source: str = "traces"            # "traces" (breaking traces vs time) | "holds"
    lowpass_hz: float = 200.0         # Butterworth cut-off
    lowpass_order: int = 2
    sg_side_points: int = 20          # Savitzky-Golay window = 2 x side points + 1
    sg_order: int = 20
    dt_s: float = 0.02                # dt of the derivative (Rashid et al. 2025: 20 ms)
    min_width: int = 50               # samples (see peak_method)
    peak_method: str = "quadratic fit"    # "quadratic fit" (LabVIEW Peak Detector) | "minimum width"
    polarity: str = "rises"           # "rises" | "drops" | "both"
    flicker_threshold: float = 0.0007
    mechanical_threshold: float = 0.002
    g_min: float = -6.5               # events outside this log G range are false positives
    g_max: float = -1.0
    level_ms: float = 20.0            # conductance before/after an event: median over this time
    t_bins: int = 200                 # conductance-time histogram of the holds
    g_bins_per_decade: int = 30

    def to_dict(self) -> dict:
        return asdict(self)


@dataclass
class Event:
    trace: int            # position in the analysed list
    index: int            # sample
    kind: str             # "flicker" | "mechanical"
    height: float         # signed dG/dt at the peak
    logG_before: float
    logG_after: float


@dataclass
class EventResult:
    events: list
    kinds: np.ndarray            # per trace: "quiet" | "flickering" | "mechanical"
    n_flicker: np.ndarray
    n_mechanical: np.ndarray
    rejected: int                # peaks dropped as false positives
    settings: EventSettings
    index: np.ndarray = field(default_factory=lambda: np.zeros(0, int))   # trace numbers

    def summary_lines(self) -> list[str]:
        n = len(self.kinds)
        if not n:
            return ["No traces."]
        c = {k: int(np.sum(self.kinds == k)) for k in ("quiet", "flickering", "mechanical")}
        s = self.settings
        return [f"{n} traces: {c['mechanical']} mechanical ({100 * c['mechanical'] / n:.1f} %), "
                f"{c['flickering']} flickering ({100 * c['flickering'] / n:.1f} %), {c['quiet']} quiet",
                f"{int(self.n_flicker.sum())} flickers, {int(self.n_mechanical.sum())} mechanical events, "
                f"{self.rejected} peaks outside log G [{s.g_min:g}, {s.g_max:g}] dropped",
                f"flickers per trace: mean {self.n_flicker.mean():.2f}"]


def savgol_smooth(y: np.ndarray, side: int, order: int) -> np.ndarray:
    """Savitzky-Golay smoothing that stays exact at high polynomial order.

    scipy's savgol_coeffs loses all precision around order 20 (the weights
    come out ~1e-17), so the least-squares projection is built here from a
    QR decomposition in the Legendre basis. Ends are mirrored.
    """
    side, order = int(side), int(order)
    w = 2 * side + 1
    if order >= w or len(y) <= w:
        return np.asarray(y, float)
    x = np.linspace(-1, 1, w)
    A = np.polynomial.legendre.legvander(x, order)
    Q, _ = np.linalg.qr(A)
    h = Q @ Q[side]                      # weights of the fitted value at the centre
    yp = np.pad(np.asarray(y, float), side, mode="reflect")
    return np.convolve(yp, h[::-1], mode="valid")


def filtered(G: np.ndarray, fs: float, s: EventSettings) -> np.ndarray:
    """Butterworth low-pass followed by Savitzky-Golay smoothing of G (in G0)."""
    y = np.asarray(G, float)
    if s.lowpass_hz and s.lowpass_hz < fs / 2 and len(y) > 3 * (2 * s.lowpass_order + 1):
        sos = butter(int(s.lowpass_order), s.lowpass_hz, fs=fs, output="sos")
        y = sosfiltfilt(sos, y)
    return savgol_smooth(y, s.sg_side_points, s.sg_order)


def derivative(y: np.ndarray, dt: float) -> np.ndarray:
    """Second-order central difference (one-sided at the ends)."""
    if len(y) < 3:
        return np.zeros(len(y))
    return np.gradient(y, dt, edge_order=2)


def find_events(d: np.ndarray, s: EventSettings):
    """Peak positions and heights (absolute) in the derivative d, per s.peak_method and s.polarity."""
    sig = {"rises": [1], "drops": [-1], "both": [1, -1]}.get(s.polarity, [1])
    pos, hts = [], []
    for sg in sig:
        y = sg * d
        if s.peak_method == "minimum width":
            pk, prop = find_peaks(y, height=s.flicker_threshold, width=s.min_width)
            pos += list(pk)
            hts += list(prop["peak_heights"])
            continue
        w = max(3, int(s.min_width) | 1)
        q = savgol_smooth(y, w // 2, 2)      # local quadratic fit over w points
        pk, prop = find_peaks(q, height=s.flicker_threshold, distance=max(1, w // 2))
        pos += list(pk)
        hts += list(prop["peak_heights"])
    order = np.argsort(pos)
    return np.asarray(pos, int)[order], np.asarray(hts, float)[order]


def detect(trace: Trace, s: EventSettings, k: int = 0):
    """Events in one trace: (events, filtered G, derivative, number of rejected peaks)."""
    G = 10.0 ** np.asarray(trace.logG, float)
    yf = filtered(G, trace.fs, s)
    d = derivative(yf, s.dt_s)
    pk, hts = find_events(d, s)
    m = max(1, int(round(s.level_ms * 1e-3 * trace.fs)))
    lg = np.log10(np.clip(yf, 1e-12, None))
    out, rejected = [], 0
    for p, h in zip(pk, hts):
        before = float(np.median(lg[max(0, p - 2 * m):max(1, p - m)]))
        after = float(np.median(lg[min(len(lg) - 1, p + m):min(len(lg), p + 2 * m)]))
        if not all(s.g_min <= v <= s.g_max for v in (before, after) if np.isfinite(v)):
            rejected += 1
            continue
        kind = "mechanical" if h >= s.mechanical_threshold else "flicker"
        out.append(Event(k, int(p), kind, float(d[p]), before, after))
    return out, yf, d, rejected


def analyse(traces: list[Trace], s: EventSettings) -> EventResult:
    events, kinds, nf, nm, rej = [], [], [], [], 0
    for k, t in enumerate(traces):
        ev, _, _, r = detect(t, s, k)
        rej += r
        events += ev
        f = sum(e.kind == "flicker" for e in ev)
        mch = sum(e.kind == "mechanical" for e in ev)
        nf.append(f)
        nm.append(mch)
        kinds.append("mechanical" if mch else "flickering" if f else "quiet")
    return EventResult(events, np.array(kinds), np.array(nf), np.array(nm), rej, s,
                       np.array([t.index for t in traces]))


# Conductance-time histogram of the holds

def _gauss(x, a, m, sg):
    return a * np.exp(-(x - m) ** 2 / (2 * sg ** 2))


def time_histogram(segments: list[Trace], s: EventSettings, g_min: float = -7.0, g_max: float = 0.0):
    """2D histogram of log G against time over all segments (counts per segment),
    plus the most probable log G and its sigma at every time bin (Gaussian fit)."""
    if not segments:
        return None
    t_max = max(len(x.logG) / x.fs for x in segments)
    te = np.linspace(0, t_max, int(s.t_bins) + 1)
    ge = np.linspace(g_min, g_max, max(4, int(round((g_max - g_min) * s.g_bins_per_decade))) + 1)
    H = np.zeros((len(te) - 1, len(ge) - 1))
    for x in segments:
        H += np.histogram2d(x.t, x.logG, [te, ge])[0]
    H /= len(segments)
    gc = 0.5 * (ge[1:] + ge[:-1])
    mp, sd = np.full(len(te) - 1, np.nan), np.full(len(te) - 1, np.nan)
    for j in range(len(te) - 1):
        h = H[j]
        if h.sum() <= 0:
            continue
        i = int(np.argmax(h))
        mp[j] = gc[i]
        try:
            popt, _ = curve_fit(_gauss, gc, h, p0=[h[i], gc[i], 0.2], maxfev=2000)
            if g_min <= popt[1] <= g_max:
                mp[j], sd[j] = popt[1], abs(popt[2])
        except (RuntimeError, ValueError):
            pass
    return 0.5 * (te[1:] + te[:-1]), gc, H, mp, sd
