"""Models for I-V curves: single-level fit and transition voltage spectroscopy.

Single-level model (SLM): transport through one molecular level at
energy eps0 from the Fermi level, broadened by the coupling Gamma to each
electrode (Gamma_L = Gamma_R = Gamma). At zero temperature the Landauer
current is

    I(V) = N G0 Gamma [atan((V/2 - e(V)) / Gamma) + atan((V/2 + e(V)) / Gamma)]

with e(V) = eps0 + a V. a is the voltage-division asymmetry (0 = the
level stays put, the bias drops equally at both contacts). Energies are
in eV and V in volts, so G0 times an energy in eV is a current in A.
The low-bias conductance is G/G0 = N Gamma^2 / (eps0^2 + Gamma^2).
Used this way in e.g. Zotti et al., Small 6, 1529 (2010).

Transition voltage spectroscopy (TVS): V_t is where V^2/|I| has its
maximum, i.e. the minimum of the Fowler-Nordheim plot ln(|I|/V^2) vs 1/V
(Beebe et al., PRL 97, 026801 (2006)). In the SLM with Gamma << eps0 the
two polarities give eps0 and a in closed form (Baldea, PRB 85, 035442
(2012)):

    eps0 = 2 |V+ V-| / sqrt(V+^2 + V-^2 + 10/3 |V+ V-|)
    a    = (V+ + V-) eps0 / (4 |V+ V-|)

and for a symmetric junction eps0 = (sqrt(3)/2) V_t.
"""

from __future__ import annotations

from dataclasses import dataclass, asdict

import numpy as np
from scipy.optimize import least_squares

from . import G0


@dataclass
class IVModelSettings:
    fit_v_max: float | None = None     # fit only |V| <= this (V); empty = the whole curve
    fit_asymmetry: bool = False        # also fit a (otherwise a = 0)
    n_molecules: float = 1.0           # N in the SLM, kept fixed
    tvs_v_min: float = 0.05            # ignore |V| below this when looking for V_t
    max_curves: int = 500              # fit at most this many single curves (all go into the mean)

    def to_dict(self) -> dict:
        return asdict(self)


def slm_current(V, eps0, Gamma, a=0.0, n=1.0):
    """SLM current in A for bias V (V), eps0 and Gamma in eV."""
    V = np.asarray(V, float)
    e = eps0 + a * V
    return n * G0 * Gamma * (np.arctan((V / 2 - e) / Gamma) + np.arctan((V / 2 + e) / Gamma))


def slm_conductance(eps0, Gamma, n=1.0):
    """Zero-bias G/G0 of the SLM."""
    return n * Gamma ** 2 / (eps0 ** 2 + Gamma ** 2)


@dataclass
class SLMFit:
    eps0: float = np.nan
    Gamma: float = np.nan
    a: float = 0.0
    eps0_err: float = np.nan
    Gamma_err: float = np.nan
    a_err: float = np.nan
    r2: float = np.nan
    ok: bool = False

    @property
    def G(self) -> float:
        return slm_conductance(self.eps0, self.Gamma) if self.ok else np.nan


def low_bias_conductance(V, I, frac: float = 0.15) -> float:
    """Slope of I(V) through |V| < frac max|V|, in G0."""
    V, I = np.asarray(V, float), np.asarray(I, float)
    sel = np.abs(V) <= frac * np.nanmax(np.abs(V))
    if np.count_nonzero(sel) < 3 or np.ptp(V[sel]) == 0:
        return np.nan
    return float(np.polyfit(V[sel], I[sel], 1)[0] / G0)


def fit_slm(V, I, s: IVModelSettings) -> SLMFit:
    """Least-squares SLM fit of one I-V curve (several starting points, best one kept)."""
    V, I = np.asarray(V, float), np.asarray(I, float)
    ok = np.isfinite(V) & np.isfinite(I)
    if s.fit_v_max is not None:
        ok &= np.abs(V) <= s.fit_v_max
    V, I = V[ok], I[ok]
    if len(V) < 8 or np.ptp(V) <= 0:
        return SLMFit()
    scale = np.max(np.abs(I)) or 1.0
    n = s.n_molecules
    g = low_bias_conductance(V, I) / n
    g = g if np.isfinite(g) and 0 < g < 0.9 else 1e-3

    def resid(p):
        a = p[2] if s.fit_asymmetry else 0.0
        return (slm_current(V, p[0], p[1], a, n) - I) / scale

    lo, hi = [0.0, 1e-6], [10.0, 5.0]
    if s.fit_asymmetry:
        lo.append(-0.5)
        hi.append(0.5)
    vmax = np.max(np.abs(V))
    best = None
    for e0 in (0.4 * vmax, 0.8 * vmax, 1.5 * vmax, 3.0 * vmax):
        G0g = min(max(e0 * np.sqrt(g / (1 - g)), 1.1e-6), 4.9)
        p0 = [min(e0, 9.9), G0g] + ([0.0] if s.fit_asymmetry else [])
        try:
            r = least_squares(resid, p0, bounds=(lo, hi), x_scale="jac")
        except ValueError:
            continue
        if best is None or r.cost < best.cost:
            best = r
    if best is None:
        return SLMFit()
    p = best.x
    dof = max(1, len(V) - len(p))
    try:
        cov = np.linalg.inv(best.jac.T @ best.jac) * (2 * best.cost / dof)
        err = np.sqrt(np.clip(np.diag(cov), 0, None))
    except np.linalg.LinAlgError:
        err = np.full(len(p), np.nan)
    ss_tot = np.sum(((I - I.mean()) / scale) ** 2)
    r2 = 1 - 2 * best.cost / ss_tot if ss_tot > 0 else np.nan
    return SLMFit(float(p[0]), float(p[1]), float(p[2]) if s.fit_asymmetry else 0.0,
                  float(err[0]), float(err[1]), float(err[2]) if s.fit_asymmetry else np.nan,
                  float(r2), True)


# Transition voltage spectroscopy

def _vt_one_side(V, I, v_min, n_bins=60):
    """Bias of the V^2/|I| maximum on one polarity (V > 0 here); NaN if it is at the edge."""
    sel = (V > v_min) & np.isfinite(I)
    if np.count_nonzero(sel) < 10:
        return np.nan
    v, i = V[sel], np.abs(I[sel])
    edges = np.linspace(v.min(), v.max(), min(n_bins, max(8, len(v) // 3)) + 1)
    k = np.clip(np.digitize(v, edges) - 1, 0, len(edges) - 2)
    cnt = np.bincount(k, minlength=len(edges) - 1)
    vb = np.bincount(k, v, len(edges) - 1)
    ib = np.bincount(k, i, len(edges) - 1)
    keep = cnt > 0
    vb, ib = vb[keep] / cnt[keep], ib[keep] / cnt[keep]
    keep = ib > 0
    vb, ib = vb[keep], ib[keep]
    if len(vb) < 5:
        return np.nan
    y = vb ** 2 / ib
    j = int(np.argmax(y))
    if j == 0 or j == len(y) - 1:
        return np.nan            # no maximum inside the measured range
    # parabola through the three points around the maximum
    x0, x1, x2 = vb[j - 1:j + 2]
    y0, y1, y2 = y[j - 1:j + 2]
    den = (x0 - x1) * (x0 - x2) * (x1 - x2)
    if den == 0:
        return float(vb[j])
    A = (x2 * (y1 - y0) + x1 * (y0 - y2) + x0 * (y2 - y1)) / den
    B = (x2 ** 2 * (y0 - y1) + x1 ** 2 * (y2 - y0) + x0 ** 2 * (y1 - y2)) / den
    xv = -B / (2 * A) if A < 0 else vb[j]
    return float(np.clip(xv, x0, x2))


def transition_voltages(V, I, v_min: float = 0.05) -> tuple[float, float]:
    """(V_t+, V_t-) of one curve; V_t- is negative. NaN where there is no maximum."""
    V, I = np.asarray(V, float), np.asarray(I, float)
    vp = _vt_one_side(V, I, v_min)
    vm = _vt_one_side(-V, -I, v_min)
    return vp, (-vm if np.isfinite(vm) else np.nan)


def eps0_from_vt(v_plus: float, v_minus: float) -> tuple[float, float]:
    """(eps0, a) from the transition voltages (Baldea 2012); one polarity -> symmetric formula."""
    if np.isfinite(v_plus) and np.isfinite(v_minus) and v_plus > 0 > v_minus:
        P = abs(v_plus * v_minus)
        eps = 2 * P / np.sqrt(v_plus ** 2 + v_minus ** 2 + 10 / 3 * P)
        return float(eps), float((v_plus + v_minus) * eps / (4 * P))
    vt = v_plus if np.isfinite(v_plus) else abs(v_minus)
    return float(np.sqrt(3) / 2 * vt) if np.isfinite(vt) else np.nan, np.nan


def fowler_nordheim(V, I, v_min: float = 0.05):
    """Points of the Fowler-Nordheim plot: (1/V, ln(|I|/V^2)) for each polarity."""
    V, I = np.asarray(V, float), np.asarray(I, float)
    out = {}
    for name, sel in (("positive", V > v_min), ("negative", V < -v_min)):
        sel = sel & (np.abs(I) > 0) & np.isfinite(I)
        out[name] = (1 / V[sel], np.log(np.abs(I[sel]) / V[sel] ** 2))
    return out


# All curves at once

@dataclass
class IVModelResult:
    index: np.ndarray
    direction: np.ndarray
    v_plus: np.ndarray
    v_minus: np.ndarray
    eps_tvs: np.ndarray
    a_tvs: np.ndarray
    g_low: np.ndarray
    fits: list                 # SLMFit per curve
    mean: dict                 # direction -> (V, I, SLMFit, v_plus, v_minus)
    settings: IVModelSettings

    def column(self, name: str) -> np.ndarray:
        return np.array([getattr(f, name) for f in self.fits], float)

    def summary_lines(self) -> list[str]:
        lines = []
        for d, (_, _, f, vp, vm) in self.mean.items():
            e, a = eps0_from_vt(vp, vm)
            lines.append(f"mean {d} curve: SLM eps0 = {f.eps0:.3f} +/- {f.eps0_err:.3f} eV, "
                         f"Gamma = {1e3 * f.Gamma:.3g} +/- {1e3 * f.Gamma_err:.2g} meV, a = {f.a:.3f}, "
                         f"R^2 = {f.r2:.4f}; V_t+ = {vp:.3f} V, V_t- = {vm:.3f} V, eps0(TVS) = {e:.3f} eV")
        eps, gam = self.column("eps0"), self.column("Gamma")
        ok = np.isfinite(eps)
        if ok.any():
            lines.append(f"single curves ({ok.sum()} fitted): median eps0 = {np.median(eps[ok]):.3f} eV, "
                         f"median Gamma = {1e3 * np.median(gam[ok]):.3g} meV")
        for name, v in (("V_t+", self.v_plus), ("V_t-", self.v_minus), ("eps0 (TVS)", self.eps_tvs)):
            v = v[np.isfinite(v)]
            if len(v):
                lines.append(f"{name}: median {np.median(v):.3f}, {len(v)} curves")
        return lines


def analyse(curves, s: IVModelSettings, mean_fn=None) -> IVModelResult:
    """TVS and SLM fit of every curve (up to s.max_curves) and of the mean curves."""
    n = len(curves)
    vp, vm, et, at, gl = (np.full(n, np.nan) for _ in range(5))
    fits = []
    for k, c in enumerate(curves):
        vp[k], vm[k] = transition_voltages(c.V, c.I, s.tvs_v_min)
        et[k], at[k] = eps0_from_vt(vp[k], vm[k])
        gl[k] = low_bias_conductance(c.V, c.I)
        fits.append(fit_slm(c.V, c.I, s) if k < s.max_curves else SLMFit())
    mean = {}
    if mean_fn is not None:
        for d, (vg, m, _) in mean_fn(curves).items():
            p, q = transition_voltages(vg, m, s.tvs_v_min)
            mean[d] = (vg, m, fit_slm(vg, m, s), p, q)
    return IVModelResult(np.array([c.index for c in curves]), np.array([c.direction for c in curves]),
                         vp, vm, et, at, gl, fits, mean, s)
