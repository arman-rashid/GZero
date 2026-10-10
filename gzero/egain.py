"""Ensemble (large-area) junctions: EGaIn and similar J(V) data.

Every file is one junction with several bias sweeps, usually
0 -> +V -> 0 -> -V -> 0. The analysis follows the statistics used for
EGaIn junctions (Reus et al., J. Phys. Chem. C 116, 6714 (2012); GaussFit by
the Chiechi group; the EGaIn module of XMe by the Hong group):

- J = I / A, with A from the contact diameter or given directly.
- A junction is a short if |J| ever reaches the short threshold; the yield is
  the fraction of junctions that are not shorts. Shorts are left out.
- Each sweep is put on a common bias grid. At every bias the values of
  log10|J| of all sweeps are histogrammed and fitted with a Gaussian: the
  centre <log|J|> and width sigma_log are what is reported. The confidence
  interval uses the number of junctions, not of sweeps, as the degrees of
  freedom: CI = t(1 - alpha/2, n_j - 1) sigma_log / sqrt(n_j - 1).
- Rectification R(V) = |J(+V)| / |J(-V)| per cycle (a set of sweeps covering
  both polarities), histogrammed as log10 R.
- V_trans from the Fowler-Nordheim minimum of every cycle.
- Optionally the single-level model on the Gaussian-mean J(V) with N free
  (molecules per cm^2 that carry the current).
- With several data sets of known molecular length d (or carbon count),
  log10 <|J|> = log10 J0 - beta d / ln 10 gives the tunnelling decay beta and J0
  (e.g. Simeone et al., JACS 135, 18131 (2013)).
"""

from __future__ import annotations

import os
import warnings
from dataclasses import dataclass, asdict, field

import numpy as np
from scipy import stats
from scipy.optimize import curve_fit

from . import bj_io
from . import ivmodels as ivm


@dataclass
class EGaInSettings:
    v_column: str = "0"
    current_column: str = "1"
    current_unit: str = "A"            # "A" | "mA" | "uA" | "nA" | "A/cm2" (already J)
    contact_diameter_um: float = 25.0  # EGaIn tip contact diameter
    contact_area_cm2: float | None = None   # overrides the diameter if set
    short_threshold: float = 1e2       # A/cm^2; |J| at or above -> the junction is a short
    v_step: float = 0.05               # bias grid (V)
    bins_per_decade: int = 10          # log|J| histograms
    alpha: float = 0.05                # confidence interval: 1 - alpha
    v_report: float = 0.5              # bias for the histograms, R and beta
    fit_slm: bool = True
    temperature: float = 300.0         # K, for the single-level fit

    def to_dict(self) -> dict:
        return asdict(self)

    def area(self) -> float:
        if self.contact_area_cm2:
            return float(self.contact_area_cm2)
        r_cm = self.contact_diameter_um * 1e-4 / 2
        return float(np.pi * r_cm ** 2)


UNIT_SCALE = {"A": 1.0, "mA": 1e-3, "uA": 1e-6, "nA": 1e-9, "A/cm2": None}


@dataclass
class Junction:
    name: str
    V: np.ndarray
    J: np.ndarray      # A/cm^2


@dataclass
class Dataset:
    name: str
    junctions: list = field(default_factory=list)
    length: float | None = None      # molecular length (nm, A or carbon number) for beta


def load_junction(path: str, s: EGaInSettings) -> Junction:
    tab = bj_io.read_table(path)
    V = tab.get(s.v_column)
    I = tab.get(s.current_column)
    if V is None or I is None:
        raise ValueError(f"{os.path.basename(path)}: bias or current column missing.")
    n = min(len(V), len(I))
    V, I = np.asarray(V[:n], float), np.asarray(I[:n], float)
    scale = UNIT_SCALE.get(s.current_unit, 1.0)
    J = I if scale is None else I * scale / s.area()
    ok = np.isfinite(V) & np.isfinite(J)
    return Junction(os.path.splitext(os.path.basename(path))[0], V[ok], J[ok])


def split_sweeps(V: np.ndarray, min_points: int = 4) -> list[tuple[int, int]]:
    """(start, stop) of every monotonic piece of the bias (EGaIn bias is a clean staircase)."""
    d = np.sign(np.diff(V))
    # carry the last non-zero direction through repeated bias values
    for i in range(1, len(d)):
        if d[i] == 0:
            d[i] = d[i - 1]
    turns = np.flatnonzero(d[1:] != d[:-1]) + 1
    edges = np.r_[0, turns, len(V) - 1]
    return [(a, b + 1) for a, b in zip(edges[:-1], edges[1:]) if b + 1 - a >= min_points]


@dataclass
class EGaInResult:
    name: str
    grid: np.ndarray                 # bias grid
    logJ: np.ndarray                 # sweeps x grid, NaN where a sweep doesn't reach
    sweep_junction: np.ndarray       # junction index of every sweep
    n_junctions: int
    n_shorts: int
    mu: np.ndarray                   # Gaussian centre of log|J| per bias
    sigma: np.ndarray
    ci: np.ndarray                   # half width of the confidence interval
    counts: np.ndarray               # sweeps per bias
    logR: dict                       # |V| -> per-cycle log10 R
    vt_plus: np.ndarray
    vt_minus: np.ndarray
    slm: ivm.SLMFit | None
    settings: EGaInSettings
    length: float | None = None

    @property
    def n_working(self) -> int:
        return self.n_junctions - self.n_shorts

    def at(self, v: float) -> int:
        return int(np.argmin(np.abs(self.grid - v)))

    def summary_lines(self) -> list[str]:
        s = self.settings
        n = self.n_junctions
        k = self.at(s.v_report)
        lines = [f"{self.name}: {n} junctions, {self.n_shorts} shorts, yield "
                 f"{100 * self.n_working / n if n else np.nan:.0f} %, {len(self.logJ)} sweeps",
                 f"at {self.grid[k]:+.2f} V: <log|J|> = {self.mu[k]:.2f}, sigma_log = {self.sigma[k]:.2f}, "
                 f"{100 * (1 - s.alpha):.0f} % CI +/- {self.ci[k]:.2f} (J in A/cm2)"]
        r = self.logR.get(round(abs(s.v_report), 6))
        if r is not None and len(r):
            mu, sg, _ = gauss_fit(r, s.bins_per_decade)
            lines.append(f"log R at {abs(s.v_report):.2f} V: {mu:.2f} +/- {sg:.2f} (R = {10 ** mu:.3g}, "
                         f"{len(r)} cycles)")
        for name, v in (("V_trans+", self.vt_plus), ("V_trans-", self.vt_minus)):
            v = v[np.isfinite(v)]
            if len(v):
                lines.append(f"{name}: median {np.median(v):.2f} V ({len(v)} cycles)")
        f = self.slm
        if f is not None and f.ok:
            lines.append(f"single-level fit of <log|J|>: eps0 = {f.eps0:.3f} +/- {f.eps0_err:.3f} eV, "
                         f"Gamma = {1e3 * f.Gamma:.3g} meV, N = {f.n:.3g} per cm2, R^2 = {f.r2:.4f}")
        return lines


def gauss_fit(values, bins_per_decade: int = 10):
    """(centre, sigma, (bin centres, counts, fit)) of a Gaussian fitted to the histogram of values."""
    v = np.asarray(values, float)
    v = v[np.isfinite(v)]
    if len(v) == 0:
        return np.nan, np.nan, None
    if len(v) < 3 or np.ptp(v) == 0:
        return float(np.mean(v)), float(np.std(v)), None
    w = 1.0 / bins_per_decade
    edges = np.arange(np.floor(v.min() / w) * w - w, v.max() + 2 * w, w)
    h, _ = np.histogram(v, edges)
    c = 0.5 * (edges[1:] + edges[:-1])
    mu, sg = float(np.mean(v)), float(np.std(v)) or w
    try:
        popt, _ = curve_fit(lambda x, a, m, s_: a * np.exp(-(x - m) ** 2 / (2 * s_ ** 2)), c, h,
                            p0=[h.max(), mu, sg], maxfev=5000)
        if v.min() - 1 <= popt[1] <= v.max() + 1:
            mu, sg = float(popt[1]), float(abs(popt[2]))
    except (RuntimeError, ValueError):
        pass
    fit = len(v) * w / (sg * np.sqrt(2 * np.pi)) * np.exp(-(c - mu) ** 2 / (2 * sg ** 2))
    return mu, sg, (c, h, fit)


def _cycles(sweeps: list[tuple[np.ndarray, np.ndarray]]):
    """Group consecutive sweeps until both polarities are covered."""
    out, cur, pos, neg = [], [], False, False
    for V, row in sweeps:
        cur.append(row)
        pos |= np.nanmax(V) > 0
        neg |= np.nanmin(V) < 0
        if pos and neg:
            with np.errstate(all="ignore"), warnings.catch_warnings():
                warnings.simplefilter("ignore", RuntimeWarning)
                out.append(np.nanmean(np.vstack(cur), axis=0) if len(cur) > 1 else cur[0])
            cur, pos, neg = [], False, False
    return out


def analyse(ds: Dataset, s: EGaInSettings) -> EGaInResult:
    if not ds.junctions:
        raise ValueError(f"Data set {ds.name} has no junctions.")
    vmax = max(np.max(np.abs(j.V)) for j in ds.junctions)
    n_grid = int(np.floor(vmax / s.v_step + 1e-9))
    grid = np.round(np.arange(-n_grid, n_grid + 1) * s.v_step, 10)
    rows, owner, cyc_rows, shorts = [], [], [], 0
    for ji, j in enumerate(ds.junctions):
        if np.any(np.abs(j.J) >= s.short_threshold):
            shorts += 1
            continue
        sw = []
        for a, b in split_sweeps(j.V):
            v, lj = j.V[a:b], np.log10(np.clip(np.abs(j.J[a:b]), 1e-300, None))
            order = np.argsort(v)
            v, lj = v[order], lj[order]
            row = np.interp(grid, v, lj, left=np.nan, right=np.nan)
            row[(grid < v.min() - 1e-9) | (grid > v.max() + 1e-9)] = np.nan
            row[np.isclose(grid, 0.0)] = np.nan          # log|J| at 0 V is meaningless
            rows.append(row)
            owner.append(ji)
            sw.append((v, row))
        cyc_rows += _cycles(sw)
    if not rows:
        raise ValueError(f"Data set {ds.name}: every junction is a short (|J| >= {s.short_threshold:g} A/cm2).")
    L = np.vstack(rows)
    owner = np.array(owner)
    n_working = len(ds.junctions) - shorts
    mu, sg, ci, cnt = (np.full(len(grid), np.nan) for _ in range(4))
    for k in range(len(grid)):
        col = L[:, k]
        col = col[np.isfinite(col)]
        cnt[k] = len(col)
        if len(col) < 2:
            continue
        mu[k], sg[k], _ = gauss_fit(col, s.bins_per_decade)
        dof = max(1, n_working - 1)
        ci[k] = stats.t.ppf(1 - s.alpha / 2, dof) * sg[k] / np.sqrt(dof)
    # rectification per cycle
    logR = {}
    C = np.vstack(cyc_rows) if cyc_rows else np.zeros((0, len(grid)))
    for k in np.flatnonzero(grid > 0):
        kn = int(np.argmin(np.abs(grid + grid[k])))
        r = C[:, k] - C[:, kn]
        logR[round(float(grid[k]), 6)] = r[np.isfinite(r)]
    # transition voltages per cycle
    vtp, vtm = [], []
    for row in C:
        ok = np.isfinite(row)
        if ok.sum() < 8:
            continue
        J = np.sign(grid[ok]) * 10 ** row[ok]
        p, q = ivm.transition_voltages(grid[ok], J, s.v_step * 0.99)
        vtp.append(p)
        vtm.append(q)
    slm = None
    if s.fit_slm:
        ok = np.isfinite(mu)
        if ok.sum() >= 8:
            ms = ivm.IVModelSettings(fit_n=True, temperature=s.temperature)
            slm = ivm.fit_slm(grid[ok], np.sign(grid[ok]) * 10 ** mu[ok], ms)
    return EGaInResult(ds.name, grid, L, owner, len(ds.junctions), shorts, mu, sg, ci, cnt, logR,
                       np.array(vtp, float), np.array(vtm, float), slm, s, ds.length)


def fit_beta(results: list[EGaInResult], v: float):
    """log10 <|J|>(v) against length over data sets: (beta, beta_err, log10 J0, log10 J0 err, points).

    beta is in 1/(unit of length): log10 J = log10 J0 - beta d / ln 10. Each point is weighted by
    the standard error of its mean, sigma_log / sqrt(n_j - 1) with n_j the working junctions.
    """
    pts = []
    for r in results:
        if r.length is None:
            continue
        k = r.at(v)
        if np.isfinite(r.mu[k]):
            sg = r.sigma[k] if np.isfinite(r.sigma[k]) and r.sigma[k] > 0 else 1.0
            sem = sg / np.sqrt(max(1, r.n_working - 1))
            pts.append((r.length, r.mu[k], sg, r.name, sem))
    if len(pts) < 2:
        return None
    d = np.array([p[0] for p in pts], float)
    y = np.array([p[1] for p in pts], float)
    w = 1 / np.array([p[4] for p in pts], float)     # weights: standard error of each mean
    A = np.vstack([np.ones_like(d), d]).T
    W = np.diag(w ** 2)
    cov = np.linalg.inv(A.T @ W @ A)
    b0, b1 = cov @ A.T @ W @ y
    if len(pts) > 2:
        res = y - (b0 + b1 * d)
        cov = cov * max(1.0, float(res @ W @ res) / (len(pts) - 2))   # never smaller than the input errors
    err = np.sqrt(np.diag(cov))
    return {"beta": float(-b1 * np.log(10)), "beta_err": float(err[1] * np.log(10)),
            "logJ0": float(b0), "logJ0_err": float(err[0]), "points": pts, "v": v}


def simulate(n_junctions: int = 20, sweeps: int = 6, beta: float = 0.9, length: float = 12.0,
             logJ0: float = 3.6, sigma: float = 0.4, rect: float = 0.0, p_short: float = 0.1,
             v_max: float = 0.5, v_step: float = 0.05, seed: int = 0, area_cm2: float = 4.9e-6) -> list[Junction]:
    """Fake EGaIn junctions with known <log|J|>, sigma_log, beta and rectification.

    log10|J|(V) = logJ0 - beta * length / ln10 + log10(V/v_max) + junction offset;
    rect = log10 R added at positive bias. Returns J in A/cm^2 (as currents I = J A).
    """
    rng = np.random.default_rng(seed)
    up = np.round(np.arange(0, v_max + v_step / 2, v_step), 10)
    cyc = np.r_[up, up[::-1][1:], -up[1:], -up[::-1][1:]]
    out = []
    for j in range(n_junctions):
        off = sigma * rng.standard_normal()
        short = rng.random() < p_short
        V = np.tile(cyc, sweeps // 2)
        mag = logJ0 - beta * length / np.log(10) + np.log10(np.clip(np.abs(V) / v_max, 1e-6, None)) \
            * 1.2 + off + rect * (V > 0) * np.abs(V) / v_max + 0.05 * rng.standard_normal(len(V))
        J = np.sign(V) * 10 ** mag
        if short:
            J = V / 1e-3 / area_cm2       # ohmic, ~1 kOhm
        out.append(Junction(f"j{j}", V, J))
    return out
