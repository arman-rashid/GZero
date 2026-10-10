"""Histograms, peak fits, plateau length and the 2D correlation map."""

from __future__ import annotations

from dataclasses import dataclass, asdict

import numpy as np
from scipy import stats
from scipy.optimize import curve_fit

from .traces import Trace


@dataclass
class HistSettings:
    g_min: float = -6.5        # log10(G/G0)
    g_max: float = 1.0
    bins_per_decade: int = 100
    z_min: float = -0.5        # nm
    z_max: float = 2.0
    z_bins: int = 200
    normalize: str = "per trace"   # "per trace" | "counts" | "density"

    def g_edges(self) -> np.ndarray:
        n = max(2, int(round((self.g_max - self.g_min) * self.bins_per_decade)))
        return np.linspace(self.g_min, self.g_max, n + 1)

    def z_edges(self) -> np.ndarray:
        return np.linspace(self.z_min, self.z_max, int(self.z_bins) + 1)

    def to_dict(self) -> dict:
        return asdict(self)


def centers(edges: np.ndarray) -> np.ndarray:
    return 0.5 * (edges[1:] + edges[:-1])


def _norm(counts: np.ndarray, n_traces: int, mode: str, width: float = 1.0) -> np.ndarray:
    if mode == "per trace" and n_traces:
        return counts / n_traces
    if mode == "density" and counts.sum():
        return counts / counts.sum() / width
    return counts.astype(float)


def hist1d(traces: list[Trace], hs: HistSettings):
    edges = hs.g_edges()
    if not traces:
        return centers(edges), np.zeros(len(edges) - 1)
    y = np.concatenate([t.logG for t in traces])
    counts, _ = np.histogram(y, edges)
    return centers(edges), _norm(counts, len(traces), hs.normalize, edges[1] - edges[0])


def hist2d(traces: list[Trace], hs: HistSettings):
    ze, ge = hs.z_edges(), hs.g_edges()
    if not traces:
        return centers(ze), centers(ge), np.zeros((len(ze) - 1, len(ge) - 1))
    z = np.concatenate([t.z for t in traces])
    y = np.concatenate([t.logG for t in traces])
    H, _, _ = np.histogram2d(z, y, [ze, ge])
    return centers(ze), centers(ge), _norm(H, len(traces), hs.normalize)


def trace_hist(trace: Trace, edges: np.ndarray) -> np.ndarray:
    return np.histogram(trace.logG, edges)[0].astype(float)


def hist_matrix(traces: list[Trace], edges: np.ndarray) -> np.ndarray:
    """Per-trace histograms, shape (n_traces, n_bins)."""
    out = np.zeros((len(traces), len(edges) - 1))
    for i, t in enumerate(traces):
        out[i] = np.histogram(t.logG, edges)[0]
    return out


# Peak fitting

def _gauss(x, a, mu, sigma):
    return a * np.exp(-0.5 * ((x - mu) / sigma) ** 2)


def _multi_gauss(x, *p):
    n = (len(p) - 1) // 3
    y = np.full_like(x, p[-1], dtype=float)
    for k in range(n):
        y += _gauss(x, *p[3 * k:3 * k + 3])
    return y


@dataclass
class PeakFit:
    center: float          # log10(G/G0)
    sigma: float           # decades
    amplitude: float
    center_err: float = np.nan
    sigma_err: float = np.nan

    @property
    def G(self) -> float:
        return 10.0 ** self.center

    @property
    def fwhm(self) -> float:
        return 2.3548200450309493 * abs(self.sigma)


def fit_peaks(x: np.ndarray, y: np.ndarray, lo: float, hi: float, n_peaks: int = 1,
              guesses: list[float] | None = None, baseline: bool = True):
    """Fit n_peaks Gaussians (+constant) to a log-conductance histogram in [lo, hi].

    Returns (list[PeakFit], x_fit, y_fit, baseline).
    """
    sel = (x >= lo) & (x <= hi)
    xs, ys = x[sel], y[sel]
    if len(xs) < 3 * n_peaks + 2 or ys.max() <= 0:
        raise ValueError("Not enough histogram data in the fit range.")
    if guesses is None:
        from scipy.signal import find_peaks
        pk, prop = find_peaks(ys, prominence=0)
        order = np.argsort(prop["prominences"])[::-1][:n_peaks] if len(pk) else []
        guesses = sorted(xs[pk[order]]) if len(order) else []
        while len(guesses) < n_peaks:
            guesses.append(lo + (hi - lo) * (len(guesses) + 1) / (n_peaks + 1))
    p0, lb, ub = [], [], []
    width = (hi - lo) / (4 * n_peaks)
    for g in guesses[:n_peaks]:
        amp = float(np.interp(g, xs, ys))
        p0 += [max(amp, ys.max() * 0.1), g, width]
        lb += [0, lo, (xs[1] - xs[0]) / 2]
        ub += [np.inf, hi, hi - lo]
    base0 = float(np.min(ys)) if baseline else 0.0
    p0.append(base0)
    lb.append(0 if baseline else -1e-12)
    ub.append(np.inf if baseline else 1e-12)
    popt, pcov = curve_fit(_multi_gauss, xs, ys, p0=p0, bounds=(lb, ub), maxfev=20000)
    perr = np.sqrt(np.clip(np.diag(pcov), 0, None))
    peaks = [PeakFit(popt[3 * k + 1], abs(popt[3 * k + 2]), popt[3 * k],
                     perr[3 * k + 1], perr[3 * k + 2]) for k in range(n_peaks)]
    peaks.sort(key=lambda p: p.center)
    xf = np.linspace(lo, hi, 400)
    return peaks, xf, _multi_gauss(xf, *popt), popt[-1]


def dwell_time_ms(peak: PeakFit, bin_width: float, fs: float, normalize: str = "per trace",
                  n_traces: int = 0) -> float:
    """Average time per trace spent in a conductance state, from the area of its fitted peak.

    tau = A / w x 1000 / fs ms (Rashid et al., JACS 147, 830 (2025), SI 2.7), where A is the
    area under the peak of the per-trace normalised 1D histogram and w the bin width:
    A / w is the number of samples per trace inside the peak.
    """
    area = peak.amplitude * abs(peak.sigma) * np.sqrt(2 * np.pi)
    if normalize == "counts" and n_traces:
        area /= n_traces
    elif normalize not in ("per trace", "counts"):
        return np.nan
    return float(area / bin_width * 1000.0 / fs) if fs > 0 and bin_width > 0 else np.nan


# Plateau length and per-trace plateau conductance

def plateau_lengths(traces: list[Trace], lo: float, hi: float, method: str = "span") -> np.ndarray:
    """Length (nm) each trace spends in the window lo <= log G <= hi.

    span  : from first point below hi to the last point above lo
                (standard 'stretching length'; robust to short excursions).
    count : number of points inside the window x step size.
    """
    out = np.zeros(len(traces))
    for i, t in enumerate(traces):
        y, z = t.logG, t.z
        if len(z) < 2:
            continue
        if method == "count":
            dz = np.median(np.abs(np.diff(z))) if len(z) > 1 else 0.0
            out[i] = np.count_nonzero((y >= lo) & (y <= hi)) * dz
            continue
        below_hi = np.flatnonzero(y <= hi)
        above_lo = np.flatnonzero(y >= lo)
        if len(below_hi) == 0 or len(above_lo) == 0:
            continue
        a, b = below_hi[0], above_lo[-1]
        out[i] = max(0.0, z[b] - z[a]) if b > a else 0.0
    return out


def describe(values: np.ndarray) -> dict:
    v = np.asarray(values, float)
    v = v[np.isfinite(v)]
    if len(v) == 0:
        return {"n": 0}
    pos = v[v > 0]
    return {
        "n": int(len(v)),
        "mean": float(np.mean(v)),
        "median": float(np.median(v)),
        "std": float(np.std(v, ddof=1)) if len(v) > 1 else 0.0,
        "variance": float(np.var(v, ddof=1)) if len(v) > 1 else 0.0,
        "skewness": float(stats.skew(v)) if len(v) > 2 else np.nan,
        "kurtosis": float(stats.kurtosis(v)) if len(v) > 3 else np.nan,
        "geometric_mean": float(stats.gmean(pos)) if len(pos) else np.nan,
        "min": float(v.min()),
        "max": float(v.max()),
    }


def plateau_conductance(traces: list[Trace], lo: float, hi: float, bins: int = 60) -> np.ndarray:
    """Most-probable log G of each trace inside [lo, hi] (NaN if no plateau)."""
    edges = np.linspace(lo, hi, bins + 1)
    c = centers(edges)
    out = np.full(len(traces), np.nan)
    for i, t in enumerate(traces):
        h = np.histogram(t.logG, edges)[0]
        if h.max() >= 5:
            out[i] = c[np.argmax(h)]
    return out


# 2D cross-correlation histogram (Makk et al., ACS Nano 2012)

def correlation2d(traces: list[Trace], edges: np.ndarray) -> np.ndarray:
    """Pearson correlation of per-trace histogram counts between all bin pairs."""
    M = hist_matrix(traces, edges)
    if len(M) < 3:
        return np.full((M.shape[1], M.shape[1]), np.nan)
    with np.errstate(divide="ignore", invalid="ignore"):
        C = np.corrcoef(M, rowvar=False)   # empty bins have zero variance -> NaN
    return np.where(np.isfinite(C), C, np.nan)


# Tunnelling decay (distance calibration helper)

def tunnelling_decay(traces: list[Trace], lo: float, hi: float) -> dict:
    """Fit log10 G = a - (beta/ln10) z in [lo, hi] for each trace.

    Returns the median decay constant beta (1/nm, in units of the current
    distance axis) and the per-trace values. With a known beta (about 20/nm
    for Au in vacuum, i.e. 2*kappa; lower in solvents) this calibrates the
    displacement ratio of an MCBJ.
    """
    betas = []
    for t in traces:
        sel = (t.logG >= lo) & (t.logG <= hi)
        if np.count_nonzero(sel) < 10:
            continue
        slope = np.polyfit(t.z[sel], t.logG[sel], 1)[0]
        if slope < 0:
            betas.append(-slope * np.log(10.0))
    betas = np.asarray(betas)
    return {"beta_median": float(np.median(betas)) if len(betas) else np.nan,
            "betas": betas, "n": int(len(betas))}
