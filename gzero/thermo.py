"""Thermopower (Seebeck coefficient) of molecular junctions.

Input: the thermoelectric voltage of the junction and the temperature, per
sample (one file per trace, or many traces in one file). The analysis is the
usual one for STM break junctions with a heated electrode (Reddy et al.,
Science 315, 1568 (2007); Widawsky et al., Nano Lett. 12, 354 (2012)):

- dT is the temperature difference across the junction. The temperature
  column is either dT itself, or the hot-side temperature with
  dT = T - T_reference. Without a temperature column every file of a data set
  gets the dT typed in for it.
- Samples (or file means) are grouped by dT in steps of dT step. At every dT
  the thermovoltage V is histogrammed and fitted with a Gaussian; its centre
  is the most probable V at that dT.
- The most probable V against dT is fitted with a straight line, weighted
  by the error of each centre. The junction Seebeck coefficient is
      S = -dV/d(dT) + S_electrode
  (S_electrode is the Seebeck coefficient of the electrode leads; 0 gives S
  relative to them). S > 0 means hole transport (HOMO closer to the Fermi
  level), S < 0 electron transport (LUMO).
- For every sample (or file) with dT != 0 also S_i = -V_i/dT_i + S_electrode is
  histogrammed, and with a conductance column it is plotted against log G.
"""

from __future__ import annotations

import os
import warnings
from dataclasses import dataclass, asdict, field

import numpy as np
from scipy.optimize import curve_fit

from . import bj_io

V_UNITS = {"V": 1e6, "mV": 1e3, "uV": 1.0, "nV": 1e-3}       # to microvolt
G_UNITS = ["log G/G0", "G/G0"]
T_MODES = ["temperature difference dT", "hot-side temperature (dT = T - T_ref)", "no column (dT per data set)"]


@dataclass
class ThermoSettings:
    v_column: str = "0"
    voltage_unit: str = "uV"
    amplifier_gain: float = 1.0          # the voltage column is divided by this
    invert_voltage: bool = False         # flips the sign of V (wiring convention)
    t_column: str = "1"
    t_mode: str = "temperature difference dT"
    t_reference: float = 295.0           # K (or degC, same unit as the column)
    g_column: str = ""                   # optional conductance column
    g_unit: str = "log G/G0"
    g_min: float | None = None           # keep samples with log G in this range
    g_max: float | None = None
    unit: str = "samples"                # histogram "samples" or "file means"
    dT_step: float = 1.0                 # K; dT is grouped to multiples of this
    min_count: int = 20                  # fewer samples / files at a dT: left out
    v_bins: int = 80
    s_electrode: float = 0.0             # uV/K added to -dV/d(dT)
    s_bins: int = 80
    bootstrap: int = 100                 # resamplings of the files for the error of each peak

    def to_dict(self) -> dict:
        return asdict(self)


@dataclass
class ThermoFile:
    name: str
    V: np.ndarray            # microvolt
    dT: np.ndarray           # K
    logG: np.ndarray | None = None


@dataclass
class ThermoDataset:
    name: str
    files: list = field(default_factory=list)


def load_file(path: str, s: ThermoSettings, dT_fixed: float | None = None) -> ThermoFile:
    tab = bj_io.read_table(path)
    V = tab.get(s.v_column)
    if V is None:
        raise ValueError(f"{os.path.basename(path)}: no voltage column.")
    V = np.asarray(V, float) * V_UNITS.get(s.voltage_unit, 1.0) / (s.amplifier_gain or 1.0)
    if s.invert_voltage:
        V = -V
    if s.t_mode.startswith("no column"):
        if dT_fixed is None:
            raise ValueError(f"{os.path.basename(path)}: no temperature column; give dT for the data set.")
        dT = np.full(len(V), float(dT_fixed))
    else:
        T = tab.get(s.t_column)
        if T is None:
            raise ValueError(f"{os.path.basename(path)}: no temperature column.")
        n = min(len(V), len(T))
        V, T = V[:n], np.asarray(T[:n], float)
        dT = T - s.t_reference if s.t_mode.startswith("hot-side") else T
    logG = None
    if s.g_column.strip():
        G = tab.get(s.g_column)
        if G is None:
            raise ValueError(f"{os.path.basename(path)}: no conductance column '{s.g_column}'.")
        n = min(len(V), len(G))
        V, dT, G = V[:n], dT[:n], np.asarray(G[:n], float)
        logG = G if s.g_unit == "log G/G0" else np.log10(np.clip(G, 1e-12, None))
    ok = np.isfinite(V) & np.isfinite(dT)
    if logG is not None:
        ok &= np.isfinite(logG)
        logG = logG[ok]
    return ThermoFile(os.path.splitext(os.path.basename(path))[0], V[ok], dT[ok], logG)


def _gauss(x, a, m, sg):
    return a * np.exp(-(x - m) ** 2 / (2 * sg ** 2))


def gauss_peak(values, bins: int):
    """(centre, sigma, error of the centre, (bin centres, counts, fitted curve))."""
    v = np.asarray(values, float)
    v = v[np.isfinite(v)]
    if len(v) < 3 or np.ptp(v) == 0:
        m = float(np.mean(v)) if len(v) else np.nan
        return m, float(np.std(v)) if len(v) else np.nan, np.nan, None
    lo, hi = np.percentile(v, [0.5, 99.5])
    if hi <= lo:
        lo, hi = v.min(), v.max()
    h, e = np.histogram(v, np.linspace(lo, hi, int(bins) + 1))
    c = 0.5 * (e[1:] + e[:-1])
    m, sg = float(np.median(v)), float(np.std(v)) or (hi - lo) / 10
    err = sg / np.sqrt(len(v))
    try:
        with warnings.catch_warnings(), np.errstate(all="ignore"):
            warnings.simplefilter("ignore")
            popt, pcov = curve_fit(_gauss, c, h, p0=[h.max(), m, sg], maxfev=5000)
        if lo <= popt[1] <= hi:
            m, sg = float(popt[1]), float(abs(popt[2]))
            # the fit error of the centre, but never below the plain standard error
            err = max(float(np.sqrt(abs(pcov[1, 1]))), sg / np.sqrt(len(v)))
    except (RuntimeError, ValueError):
        pass
    return m, sg, err, (c, h, _gauss(c, h.max(), m, sg) if np.isfinite(sg) and sg > 0 else h * 0)


@dataclass
class ThermoResult:
    name: str
    dT_levels: np.ndarray        # grouped dT values
    v_peak: np.ndarray           # most probable V per dT (uV)
    v_sigma: np.ndarray
    v_err: np.ndarray
    counts: np.ndarray
    groups: dict                 # dT -> V values (uV)
    slope: float                 # uV/K
    slope_err: float
    intercept: float
    S: float                     # uV/K
    S_err: float
    S_values: np.ndarray         # per sample / file
    S_logG: np.ndarray | None    # log G of those, if known
    S_peak: float
    S_sigma: float
    settings: ThermoSettings
    n_files: int = 0

    def summary_lines(self) -> list[str]:
        s = self.settings
        lines = [f"{self.name}: {self.n_files} files, {len(self.dT_levels)} temperature differences, "
                 f"{int(self.counts.sum())} {s.unit}"]
        for d, m, e, n in zip(self.dT_levels, self.v_peak, self.v_err, self.counts):
            lines.append(f"dT = {d:g} K: V = {m:.3f} +/- {e:.3f} uV ({int(n)})")
        if np.isfinite(self.S):
            kind = "hole transport (HOMO)" if self.S > 0 else "electron transport (LUMO)"
            lines.append(f"<b>S = {self.S:.2f} +/- {self.S_err:.2f} uV/K</b> from the slope of V vs dT "
                         f"(dV/d(dT) = {self.slope:.3f} uV/K, offset {self.intercept:.3f} uV); {kind}")
        if np.isfinite(self.S_peak):
            lines.append(f"S histogram: peak {self.S_peak:.2f} uV/K, width {self.S_sigma:.2f} uV/K "
                         f"({len(self.S_values)} {s.unit} with dT != 0)")
        return lines


def analyse(ds: ThermoDataset, s: ThermoSettings) -> ThermoResult:
    if not ds.files:
        raise ValueError(f"Data set {ds.name} has no files.")
    V, D, Gs, F = [], [], [], []
    for fi, f in enumerate(ds.files):
        keep = np.ones(len(f.V), bool)
        if f.logG is not None and (s.g_min is not None or s.g_max is not None):
            if s.g_min is not None:
                keep &= f.logG >= s.g_min
            if s.g_max is not None:
                keep &= f.logG <= s.g_max
        if not keep.any():
            continue
        if s.unit == "file means":
            # one value per file and dT level
            lev = np.round(f.dT[keep] / s.dT_step) * s.dT_step
            for L in np.unique(lev):
                sel = lev == L
                V.append([np.mean(f.V[keep][sel])])
                D.append([L])
                Gs.append([np.median(f.logG[keep][sel])] if f.logG is not None else [np.nan])
                F.append([fi])
        else:
            V.append(f.V[keep])
            D.append(f.dT[keep])
            Gs.append(f.logG[keep] if f.logG is not None else np.full(keep.sum(), np.nan))
            F.append(np.full(keep.sum(), fi))
    if not V:
        raise ValueError(f"Data set {ds.name}: no data left (check the conductance window).")
    V, D, Gs, F = np.concatenate(V), np.concatenate(D), np.concatenate(Gs), np.concatenate(F)
    lev = np.round(D / s.dT_step) * s.dT_step
    levels, groups, owners = [], {}, {}
    for L in np.unique(lev):
        sel = lev == L
        if sel.sum() >= max(1, s.min_count):
            levels.append(float(L))
            groups[float(L)] = V[sel]
            owners[float(L)] = F[sel]
    if not levels:
        raise ValueError(f"Data set {ds.name}: fewer than {s.min_count} {s.unit} at every dT.")
    pk = [gauss_peak(groups[L], s.v_bins) for L in levels]
    vp = np.array([p[0] for p in pk])
    vs = np.array([p[1] for p in pk])
    # samples of one file are not independent, and a Gaussian fit to a histogram scatters more than
    # sigma/sqrt(n): the error of each peak comes from resampling whole files (bootstrap)
    ve = np.array([p[2] for p in pk])
    rng = np.random.default_rng(0)
    for k, L in enumerate(levels):
        ids = np.unique(owners[L])
        if s.bootstrap < 10 or len(ids) < 3:
            continue
        by_file = {i: groups[L][owners[L] == i] for i in ids}
        boots = []
        for _ in range(int(s.bootstrap)):
            pick = rng.choice(ids, len(ids), replace=True)
            boots.append(gauss_peak(np.concatenate([by_file[i] for i in pick]), s.v_bins)[0])
        boots = np.asarray(boots)
        if np.isfinite(boots).sum() >= 10:
            ve[k] = max(ve[k], float(np.nanstd(boots)))
    cnt = np.array([len(groups[L]) for L in levels])
    dl = np.array(levels)
    slope = slope_err = intercept = S = S_err = np.nan
    if len(levels) >= 2:
        w = 1 / np.where(np.isfinite(ve) & (ve > 0), ve, np.nanmax(vs) or 1.0)
        A = np.vstack([np.ones_like(dl), dl]).T
        W = np.diag(w ** 2)
        cov = np.linalg.inv(A.T @ W @ A)
        b0, b1 = cov @ A.T @ W @ vp
        if len(levels) > 2:
            res = vp - (b0 + b1 * dl)
            cov = cov * max(1.0, float(res @ W @ res) / (len(levels) - 2))
        intercept, slope = float(b0), float(b1)
        slope_err = float(np.sqrt(cov[1, 1]))
        S, S_err = -slope + s.s_electrode, slope_err
    elif levels[0] != 0:
        S, S_err = -vp[0] / dl[0] + s.s_electrode, ve[0] / abs(dl[0])
    nz = np.abs(D) >= s.dT_step / 2
    Si = -V[nz] / D[nz] + s.s_electrode
    Sg = Gs[nz] if np.isfinite(Gs).any() else None
    sp, ss, _, _ = gauss_peak(Si, s.s_bins) if len(Si) >= 3 else (np.nan, np.nan, np.nan, None)
    return ThermoResult(ds.name, dl, vp, vs, ve, cnt, groups, slope, slope_err, intercept, S, S_err,
                        Si, Sg, sp, ss, s, len(ds.files))


def simulate(S: float = 10.0, dTs=(0.0, 10.0, 20.0, 30.0), traces: int = 40, points: int = 200,
             noise_uV: float = 1.5, spread: float = 3.0, s_electrode: float = 0.0, seed: int = 0,
             logG: float = -3.0) -> list[ThermoFile]:
    """Fake thermovoltage traces with a known junction Seebeck coefficient S (uV/K).

    Every trace has its own S drawn with the given spread (junction-to-junction variation)
    and white noise on V; log G scatters by 0.2 decades around logG."""
    rng = np.random.default_rng(seed)
    out = []
    for dT in dTs:
        for k in range(traces):
            s_k = S + spread * rng.standard_normal()
            V = -(s_k - s_electrode) * dT + noise_uV * rng.standard_normal(points)
            out.append(ThermoFile(f"dT{dT:g}_{k}", V, np.full(points, float(dT)),
                                  logG + 0.2 * rng.standard_normal(points)))
    return out
