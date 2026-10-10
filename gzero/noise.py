"""Flicker noise.

NP = PSD of G integrated from f_lo to f_hi (G in G0, so NP is in G0^2).
NP ~ G^n with n ~ 1 for through-bond and n ~ 2 for through-space coupling
(Adak 2015).

Two ways to get n:
- the bell: 2D histogram of log(NP/G^n) vs log G, fitted with a 2D
  Gaussian; n is where the fit has no correlation (Adak 2015)
- Morris et al., JPCC 2025: trim 5 ms at both ends, only keep traces whose
  ends sit on the molecular peak, ADF test on log|G|, subtract the mean
  before the PSD, and take the Theil-Sen slope instead of least squares
"""

from __future__ import annotations

import warnings
from dataclasses import dataclass, asdict, field

import numpy as np
from scipy import signal, stats
from scipy.integrate import trapezoid
from scipy.optimize import curve_fit

from .stationarity import adfuller
from .traces import Trace


@dataclass
class NoiseSettings:
    source: str = "holds"          # "holds" (piezo-hold periods) | "plateau windows" (in traces)
    cut_initial_ms: float = 10.0   # skip mechanical relaxation after the piezo stops
    cut_final_ms: float = 0.0      # skip the end of the hold
    window_ms: float = 100.0       # analysis window; 0 = whole (trimmed) hold
    max_windows_per_segment: int = 0   # 0 = all
    f_lo: float = 100.0            # integrate PSD from (Hz)
    f_hi: float = 1000.0           # ... to (Hz)
    psd_method: str = "welch"      # "welch" | "periodogram"
    psd_window: str = "hann"       # "hann" | "boxcar" (rectangular, plain one-sided spectrum)
    welch_nperseg: int = 256
    detrend: str = "constant"      # "constant" (subtract mean) | "linear"
    # selection
    g_min: float = -5.0            # log10(G/G0) of window mean
    g_max: float = -1.0
    max_abs_kurtosis: float | None = None    # excess kurtosis limit (rejects switching); None = off
    max_drift: float | None = 0.3            # |log G(first) - log G(last)| limit (decades)
    drift_points: int = 50                   # averaging points at each end for the drift
    stationarity_test: bool = False          # Augmented Dickey-Fuller on log10|G|
    adf_alpha: float = 0.05                  # keep if p <= alpha (stationary)
    peak_center: float | None = None         # log G of the molecular peak (None = off)
    peak_sigma: float | None = None          # its Gaussian sigma (decades)
    peak_nsigma: float = 2.0                 # ends must lie within centre +/- nsigma * sigma
    end_fraction: float = 0.05               # "ends" = first and last fraction of the window
    end_points: int = 0                      # if > 0: "ends" = this many points instead of the fraction
    # instrument floor
    subtract_floor: bool = False
    floor_below: float = -5.2      # windows with log G below this define the floor PSD
    # scaling exponent
    estimator: str = "Theil-Sen"   # "Theil-Sen" | "OLS" | "2D Gaussian"
    n_min: float = 0.5
    n_max: float = 3.0
    n_step: float = 0.05
    x_bins: int = 40
    y_bins: int = 40
    bootstrap: int = 300

    def to_dict(self) -> dict:
        return asdict(self)


def morris2025_settings(base: NoiseSettings | None = None) -> NoiseSettings:
    """Settings of Morris et al., J. Phys. Chem. C 2025 (keeps the user's ranges/peak)."""
    s = NoiseSettings(**asdict(base)) if base is not None else NoiseSettings()
    s.cut_initial_ms = 5.0
    s.cut_final_ms = 5.0
    s.window_ms = 0.0              # one noise value per (trimmed) hold
    s.f_lo, s.f_hi = 100.0, 1000.0
    s.psd_method = "periodogram"
    s.psd_window = "boxcar"
    s.detrend = "constant"
    s.max_abs_kurtosis = None
    s.max_drift = None
    s.stationarity_test = True
    s.adf_alpha = 0.05
    s.peak_nsigma = 2.0
    s.end_fraction = 0.05
    s.estimator = "Theil-Sen"
    return s


def rashid2025_settings(base: NoiseSettings | None = None) -> NoiseSettings:
    """Protocol of Rashid et al., JACS 147, 830 (2025), SI 2.3 (and JACS 146, 9063 (2024)).

    Piezo held 160 ms at 100 mV; the first 10 ms are dropped and the rest is one
    window. The first and last 100 points must average within +/- 1 sigma of the most
    probable conductance. DFT squared (periodogram), integrated 100 Hz - 1 kHz.
    n is where the Pearson r of NP/G^n and G is smallest, scanned 0.3 - 2.3 in 0.01.
    """
    s = NoiseSettings(**asdict(base)) if base is not None else NoiseSettings()
    s.cut_initial_ms, s.cut_final_ms = 10.0, 0.0
    s.window_ms = 0.0
    s.f_lo, s.f_hi = 100.0, 1000.0
    s.psd_method, s.psd_window = "periodogram", "boxcar"
    s.detrend = "constant"
    s.max_abs_kurtosis = None
    s.max_drift = None
    s.stationarity_test = False
    s.peak_nsigma = 1.0
    s.end_points = 100
    s.estimator = "OLS"
    s.n_min, s.n_max, s.n_step = 0.3, 2.3, 0.01
    return s


@dataclass
class NoiseWindow:
    segment: int            # index of the source hold / trace
    start: int              # sample offset inside the segment
    length: int             # samples
    G: float                # mean G/G0
    NP: float               # integrated noise power (G0^2)
    kurtosis: float
    drift: float
    adf_p: float = np.nan
    accepted: bool = True
    reason: str = ""


@dataclass
class NoiseResult:
    windows: list
    freqs: np.ndarray
    psd_mean: np.ndarray            # mean PSD of accepted windows
    psd_floor: np.ndarray | None
    alpha: float = np.nan           # PSD ~ 1/f^alpha in [f_lo, f_hi]
    settings: dict = field(default_factory=dict)

    def accepted(self):
        return [w for w in self.windows if w.accepted]

    def arrays(self):
        acc = self.accepted()
        return np.array([w.G for w in acc]), np.array([w.NP for w in acc])

    def rejected_counts(self) -> dict:
        out = {}
        for w in self.windows:
            if not w.accepted:
                out[w.reason] = out.get(w.reason, 0) + 1
        return out


def _psd(x: np.ndarray, fs: float, ns: NoiseSettings):
    if ns.psd_method == "periodogram":
        return signal.periodogram(x, fs, detrend=ns.detrend, window=ns.psd_window)
    nper = min(len(x), ns.welch_nperseg)
    return signal.welch(x, fs, nperseg=nper, detrend=ns.detrend, window=ns.psd_window)


def _windows_from(segment: Trace, ns: NoiseSettings):
    fs = segment.fs
    y = segment.logG
    a0 = int(round(ns.cut_initial_ms * 1e-3 * fs))
    a1 = len(y) - int(round(ns.cut_final_ms * 1e-3 * fs))
    y = y[a0:a1]
    if len(y) < 16:
        return []
    if ns.window_ms <= 0:
        return [(a0, y)]
    w = int(round(ns.window_ms * 1e-3 * fs))
    if w < 16 or len(y) < w:
        return []
    out = []
    for k in range(len(y) // w):
        out.append((a0 + k * w, y[k * w:(k + 1) * w]))
        if ns.max_windows_per_segment and len(out) >= ns.max_windows_per_segment:
            break
    return out


def _check(win: NoiseWindow, y: np.ndarray, ns: NoiseSettings):
    """Apply the selection criteria in order; set accepted/reason (and adf_p)."""
    ylog = np.log10(win.G)
    if not (ns.g_min <= ylog <= ns.g_max):
        win.accepted, win.reason = False, "conductance out of range"
        return
    if ns.peak_center is not None and ns.peak_sigma is not None:
        m = int(ns.end_points) if ns.end_points > 0 else max(1, int(round(ns.end_fraction * len(y))))
        m = min(m, max(1, len(y) // 2))
        lo = ns.peak_center - ns.peak_nsigma * ns.peak_sigma
        hi = ns.peak_center + ns.peak_nsigma * ns.peak_sigma
        g = 10.0 ** y
        ends = (np.log10(np.mean(g[:m])), np.log10(np.mean(g[-m:])))
        if not all(lo <= e <= hi for e in ends):
            win.accepted, win.reason = False, "ends outside peak"
            return
    if ns.max_abs_kurtosis is not None and abs(win.kurtosis) > ns.max_abs_kurtosis:
        win.accepted, win.reason = False, "kurtosis"
        return
    if ns.max_drift is not None and win.drift > ns.max_drift:
        win.accepted, win.reason = False, "drift"
        return
    if ns.stationarity_test:
        try:
            win.adf_p = adfuller(y).pvalue          # y = log10|G|
        except ValueError:
            win.adf_p = np.nan
        if not (win.adf_p <= ns.adf_alpha):
            win.accepted, win.reason = False, "non-stationary (ADF)"


def analyse(segments: list[Trace], ns: NoiseSettings) -> NoiseResult:
    """Compute noise power for every window of every segment and apply selection."""
    windows, psds, floor_psds = [], [], []
    freqs = None
    for si, seg in enumerate(segments):
        if ns.source == "plateau windows":
            inside = (seg.logG >= ns.g_min) & (seg.logG <= ns.g_max)
            if np.count_nonzero(inside) < 16:
                continue
        for start, y in _windows_from(seg, ns):
            g = 10.0 ** y
            f, p = _psd(g, seg.fs, ns)
            if freqs is None:
                freqs = f
            elif len(f) != len(freqs) or not np.allclose(f, freqs):
                # windows of different length (whole holds): put on a common grid
                p = np.interp(freqs, f, p)
                f = freqs
            band = (f >= ns.f_lo) & (f <= ns.f_hi)
            if np.count_nonzero(band) < 2:
                raise ValueError("The integration band contains fewer than 2 frequency points. Widen it or use longer windows.")
            NP = float(trapezoid(p[band], f[band]))
            gm = float(np.mean(g))
            with np.errstate(all="ignore"), warnings.catch_warnings():
                warnings.simplefilter("ignore")   # flat windows: kurtosis undefined
                k = float(stats.kurtosis(g)) if len(g) > 3 else np.nan
            n = min(ns.drift_points, len(y) // 4) or 1
            drift = float(abs(np.mean(y[:n]) - np.mean(y[-n:])))
            win = NoiseWindow(si, start, len(y), gm, NP, k, drift)
            if np.log10(gm) < ns.floor_below:
                floor_psds.append(p)
            _check(win, y, ns)
            windows.append(win)
            psds.append(p if win.accepted else None)

    if freqs is None:
        return NoiseResult([], np.array([]), np.array([]), None, settings=ns.to_dict())

    floor = np.mean(floor_psds, axis=0) if floor_psds else None
    if ns.subtract_floor and floor is not None:
        band = (freqs >= ns.f_lo) & (freqs <= ns.f_hi)
        floor_np = float(trapezoid(floor[band], freqs[band]))
        for w in windows:
            w.NP = max(w.NP - floor_np, np.finfo(float).tiny)

    acc = [p for p in psds if p is not None]
    psd_mean = np.mean(acc, axis=0) if acc else np.zeros_like(freqs)
    alpha = np.nan
    band = (freqs >= ns.f_lo) & (freqs <= ns.f_hi) & (psd_mean > 0)
    if np.count_nonzero(band) > 3:
        alpha = -np.polyfit(np.log10(freqs[band]), np.log10(psd_mean[band]), 1)[0]
    return NoiseResult(windows, freqs, psd_mean, floor, float(alpha), ns.to_dict())


# Scaling exponent n

def _gauss2d(xy, A, x0, y0, sx, sy, rho, c):
    x, y = xy
    dx, dy = (x - x0) / sx, (y - y0) / sy
    q = (dx * dx - 2 * rho * dx * dy + dy * dy) / (2 * (1 - rho * rho))
    return (A * np.exp(-q) + c).ravel()


@dataclass
class BellFit:
    n: float
    x0: float          # log10 G centre
    y0: float          # log10 (NP/G^n) centre
    sx: float
    sy: float
    rho: float         # correlation of the fitted Gaussian
    amplitude: float
    ok: bool = True


def fit_bell(x: np.ndarray, y: np.ndarray, n: float, xb: int = 40, yb: int = 40):
    """Fit a 2D Gaussian to the histogram of (x=log G, y=log NP/G^n)."""
    H, xe, ye = np.histogram2d(x, y, [xb, yb])
    xc, yc = 0.5 * (xe[1:] + xe[:-1]), 0.5 * (ye[1:] + ye[:-1])
    X, Y = np.meshgrid(xc, yc, indexing="ij")
    sx0, sy0 = np.std(x) or 0.1, np.std(y) or 0.1
    rho0 = np.clip(np.corrcoef(x, y)[0, 1], -0.9, 0.9) if len(x) > 2 else 0.0
    p0 = [H.max(), np.mean(x), np.mean(y), sx0, sy0, rho0, 0.0]
    lb = [0, xe[0], ye[0], 1e-3, 1e-3, -0.99, 0]
    ub = [np.inf, xe[-1], ye[-1], 10 * (xe[-1] - xe[0]), 10 * (ye[-1] - ye[0]), 0.99, H.max()]
    try:
        popt, _ = curve_fit(_gauss2d, (X, Y), H.ravel(), p0=p0, bounds=(lb, ub), maxfev=20000)
        fit = BellFit(n, popt[1], popt[2], popt[3], popt[4], popt[5], popt[0])
    except (RuntimeError, ValueError):
        fit = BellFit(n, p0[1], p0[2], p0[3], p0[4], p0[5], p0[0], ok=False)
    return fit, H, xe, ye


def theil_sen(x: np.ndarray, y: np.ndarray, max_points: int = 6000, seed: int = 0) -> dict:
    """Theil-Sen slope with Sen's 95 % confidence interval and standard error.

    More than max_points points are randomly subsampled (the estimator
    needs all N(N-1)/2 pairwise slopes).
    """
    if len(x) > max_points:
        idx = np.random.default_rng(seed).choice(len(x), max_points, replace=False)
        x, y = x[idx], y[idx]
    slope, intercept, lo, hi = stats.theilslopes(y, x, alpha=0.95)
    se = (hi - lo) / (2 * 1.959963984540054)
    return {"slope": float(slope), "intercept": float(intercept), "ci95": (float(lo), float(hi)),
            "se": float(se), "n_used": int(len(x))}


def scaling_exponent(G: np.ndarray, NP: np.ndarray, ns: NoiseSettings) -> dict:
    """Determine the noise-scaling exponent n (NP ~ G^n) with every estimator."""
    ok = (G > 0) & (NP > 0)
    x = np.log10(G[ok])
    ly = np.log10(NP[ok])
    if len(x) < 5:
        raise ValueError("At least 5 accepted noise windows are needed.")
    grid = np.arange(ns.n_min, ns.n_max + 1e-9, ns.n_step)

    # OLS = zero of the Pearson correlation of (log G, log NP/G^n)
    lr = stats.linregress(x, ly)
    n_ols, se_ols = float(lr.slope), float(lr.stderr)
    rho_p = np.array([np.corrcoef(x, ly - n * x)[0, 1] for n in grid])
    n_pearson_grid = float(grid[np.nanargmin(np.abs(rho_p))])   # Rashid et al. 2025: minimum |r| on the grid

    # Theil-Sen = zero of Kendall's tau of (log G, log NP/G^n)
    ts = theil_sen(x, ly)
    sub = slice(None) if len(x) <= 3000 else np.random.default_rng(1).choice(len(x), 3000, replace=False)
    tau_k = np.array([stats.kendalltau(x[sub], (ly - n * x)[sub])[0] for n in grid])

    # 2D Gaussian ("bell") fit: zero of the fitted correlation
    rho_fit = []
    for n in grid:
        f, *_ = fit_bell(x, ly - n * x, n, ns.x_bins, ns.y_bins)
        rho_fit.append(f.rho if f.ok else np.nan)
    rho_fit = np.array(rho_fit)
    n_fit = _zero_crossing(grid, rho_fit)

    # bootstrap check of the OLS uncertainty
    rng = np.random.default_rng(0)
    boots = []
    for _ in range(int(ns.bootstrap)):
        idx = rng.integers(0, len(x), len(x))
        vx = np.var(x[idx], ddof=1)
        if vx > 0:
            boots.append(np.cov(x[idx], ly[idx])[0, 1] / vx)
    ci_boot = tuple(np.percentile(boots, [2.5, 97.5])) if boots else (np.nan, np.nan)

    if ns.estimator == "OLS":
        best = n_ols
    elif ns.estimator == "2D Gaussian" and np.isfinite(n_fit):
        best = n_fit
    else:
        best = ts["slope"]
    bell, H, xe, ye = fit_bell(x, ly - best * x, best, ns.x_bins, ns.y_bins)
    return {
        "n_grid": grid, "rho_pearson": rho_p, "tau_kendall": tau_k, "rho_fit": rho_fit,
        "n_pearson_grid": n_pearson_grid,
        "n_ols": n_ols, "n_ols_se": se_ols, "n_ols_ci95": ci_boot, "ols_intercept": float(lr.intercept),
        "n_tse": ts["slope"], "n_tse_se": ts["se"], "n_tse_ci95": ts["ci95"], "tse_intercept": ts["intercept"],
        "tse_points": ts["n_used"],
        "n_fit": float(n_fit), "n_best": float(best), "estimator": ns.estimator,
        "bell": bell, "hist": H, "x_edges": xe, "y_edges": ye,
        "x": x, "log_np": ly, "n_points": int(len(x)),
        "interpretation": interpret(best),
    }


def _zero_crossing(xs: np.ndarray, ys: np.ndarray) -> float:
    good = np.isfinite(ys)
    xs, ys = xs[good], ys[good]
    s = np.flatnonzero(np.sign(ys[:-1]) * np.sign(ys[1:]) <= 0)
    if len(s) == 0:
        return np.nan
    i = s[0]
    if ys[i + 1] == ys[i]:
        return float(xs[i])
    return float(xs[i] - ys[i] * (xs[i + 1] - xs[i]) / (ys[i + 1] - ys[i]))


def interpret(n: float) -> str:
    if not np.isfinite(n):
        return "undetermined"
    if n < 1.3:
        return "n ~ 1: through-bond coupling dominates"
    if n < 1.7:
        return "1.3 < n < 1.7: mixed through-bond / through-space coupling"
    return "n ~ 2: through-space (direct tunnelling) coupling dominates"


def binned_noise(G: np.ndarray, NP: np.ndarray, bins: int = 20):
    """Mean log NP in log G bins (for overlaying a trend on the scatter)."""
    x, y = np.log10(G), np.log10(NP)
    e = np.linspace(x.min(), x.max(), bins + 1)
    idx = np.digitize(x, e) - 1
    xc, ym, ys = [], [], []
    for k in range(bins):
        s = idx == k
        if np.count_nonzero(s) >= 3:
            xc.append(0.5 * (e[k] + e[k + 1]))
            ym.append(np.mean(y[s]))
            ys.append(np.std(y[s]))
    return np.array(xc), np.array(ym), np.array(ys)
