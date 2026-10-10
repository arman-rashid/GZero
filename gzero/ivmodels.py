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

Finite temperature: the Fermi functions of the leads smear the bias
window. Two ways to compute it, which agree to better than 1e-9:

- analytic: with L(E) = Gamma^2 / (E^2 + Gamma^2),
  int L(E) f(E - mu) dE = pi Gamma [1/2 - Im psi(1/2 + (Gamma - i mu) / (2 pi kT)) / pi]
  (psi = digamma). Fast, used for the single curves.
- energy grid: the Landauer integral done numerically on a grid of
  energies, as in Rashid et al., JACS 147, 830 (2025), SI section 2.4
  (step 0.02 meV, -10 to 10 eV, Levenberg-Marquardt fit of the
  most probable I-V curve).
"""

from __future__ import annotations

from dataclasses import dataclass, asdict

import numpy as np
from scipy.optimize import least_squares
from scipy.signal import fftconvolve
from scipy.special import psi

from . import G0

KB_EV = 8.617333262e-5      # Boltzmann constant in eV/K


@dataclass
class IVModelSettings:
    fit_v_max: float | None = None     # fit only |V| <= this (V); empty = the whole curve
    fit_asymmetry: bool = False        # also fit a (otherwise a = 0)
    n_molecules: float = 1.0           # N in the SLM, kept fixed
    tvs_v_min: float = 0.05            # ignore |V| below this when looking for V_t
    max_curves: int = 500              # fit at most this many single curves (all go into the mean)
    temperature: float = 0.0           # K; 0 = zero-temperature closed form
    integration: str = "analytic"      # "analytic" | "energy grid" (mean / most probable curve)
    energy_step: float = 0.02          # meV, energy grid step (Rashid et al. 2025: 0.02 meV)
    energy_min: float = -10.0          # eV, energy grid limits (Rashid et al. 2025: -10 .. 10 eV)
    energy_max: float = 10.0
    fit_method: str = "trust region"   # "trust region" (bounded) | "Levenberg-Marquardt"
    fit_n: bool = False                # fit N too (large-area / ensemble junctions)
    mp_v_bins: int = 40                # most probable curve: bias bins ...
    mp_i_bins: int = 100               # ... and current bins of the 2D I-V histogram
    weight_by_sigma: bool = False      # weight the most probable curve by 1/sigma of each bias bin

    def to_dict(self) -> dict:
        return asdict(self)


def rashid2025_settings(base: IVModelSettings | None = None) -> IVModelSettings:
    """I-V fit as in Rashid et al., JACS 147, 830 (2025), SI 2.4: Landauer integral with Fermi functions
    on an energy grid (-10 .. 10 eV, 0.02 meV), Gamma_L = Gamma_R, Levenberg-Marquardt fit of the
    most probable I-V curve. The SI gives no temperature; 300 K is assumed."""
    s = IVModelSettings(**asdict(base)) if base is not None else IVModelSettings()
    s.temperature = 300.0
    s.integration = "energy grid"
    s.energy_step, s.energy_min, s.energy_max = 0.02, -10.0, 10.0
    s.fit_method = "Levenberg-Marquardt"
    s.fit_asymmetry = False
    s.fit_n = False
    s.n_molecules = 1.0
    return s


def _occupied(y, Gamma, kT):
    """int L(u) f(u - y) du with L(u) = Gamma^2/(u^2 + Gamma^2); f has chemical potential y."""
    if kT <= 0:
        return Gamma * (np.pi / 2 + np.arctan(y / Gamma))
    z = 0.5 + (Gamma - 1j * np.asarray(y, float)) / (2 * np.pi * kT)
    return np.pi * Gamma * (0.5 - np.imag(psi(z)) / np.pi)


def slm_current(V, eps0, Gamma, a=0.0, n=1.0, T=0.0):
    """SLM current in A for bias V (V), eps0 and Gamma in eV, temperature T in K.

    Symmetric bias (mu_L = V/2, mu_R = -V/2); the level sits at eps0 + a V.
    """
    V = np.asarray(V, float)
    Gamma = abs(Gamma)
    e = abs(eps0) + a * V
    kT = KB_EV * T
    return n * G0 * (_occupied(V / 2 - e, Gamma, kT) - _occupied(-V / 2 - e, Gamma, kT))


def slm_current_grid(V, eps0, Gamma, a=0.0, n=1.0, T=0.0, step_meV=0.02, e_min=-10.0, e_max=10.0):
    """The same current from the Landauer integral on an energy grid (Rashid et al. 2025, SI 2.4).

    I = (2e/h) int Tr(E) [f_L(E) - f_R(E)] dE on E = e_min .. e_max in steps of step_meV.
    Tr is integrated once (C(E) = int Tr); the Fermi functions enter by convolving C
    with -df/dE, so one evaluation costs one FFT whatever the number of bias points.
    """
    V = np.asarray(V, float)
    Gamma, eps0 = abs(Gamma), abs(eps0)
    h = step_meV * 1e-3
    E = np.arange(e_min, e_max + h / 2, h)
    tr = Gamma ** 2 / ((E - eps0) ** 2 + Gamma ** 2)
    C = np.concatenate([[0.0], np.cumsum(0.5 * (tr[1:] + tr[:-1]) * h)])
    kT = KB_EV * T
    if kT > 0:
        half = min(int(np.ceil(40 * kT / h)), len(E) // 2 - 1)
        u = np.arange(-half, half + 1) * h
        w = 1.0 / (4 * kT * np.cosh(u / (2 * kT)) ** 2) * h
        C = fftconvolve(C, w / w.sum(), mode="same")
    # level at eps0 + a V: same as shifting both chemical potentials by -a V
    mu_l, mu_r = V / 2 - a * V, -V / 2 - a * V
    return n * G0 * (np.interp(mu_l, E, C) - np.interp(mu_r, E, C))


def model_current(V, eps0, Gamma, a, n, s: "IVModelSettings", grid: bool = False):
    if grid and s.integration == "energy grid":
        return slm_current_grid(V, eps0, Gamma, a, n, s.temperature, s.energy_step, s.energy_min, s.energy_max)
    return slm_current(V, eps0, Gamma, a, n, s.temperature)


def slm_conductance(eps0, Gamma, n=1.0):
    """Zero-bias G/G0 of the SLM."""
    return n * Gamma ** 2 / (eps0 ** 2 + Gamma ** 2)


@dataclass
class SLMFit:
    eps0: float = np.nan
    Gamma: float = np.nan
    a: float = 0.0
    n: float = 1.0
    eps0_err: float = np.nan
    Gamma_err: float = np.nan
    a_err: float = np.nan
    r2: float = np.nan
    ok: bool = False
    n_err: float = np.nan

    @property
    def G(self) -> float:
        return slm_conductance(self.eps0, self.Gamma, self.n) if self.ok else np.nan


def low_bias_conductance(V, I, frac: float = 0.15) -> float:
    """Slope of I(V) through |V| < frac max|V|, in G0."""
    V, I = np.asarray(V, float), np.asarray(I, float)
    sel = np.abs(V) <= frac * np.nanmax(np.abs(V))
    if np.count_nonzero(sel) < 3 or np.ptp(V[sel]) == 0:
        return np.nan
    return float(np.polyfit(V[sel], I[sel], 1)[0] / G0)


def fit_slm(V, I, s: IVModelSettings, sigma=None, grid: bool = False) -> SLMFit:
    """SLM fit of one I-V curve (several starting points, best one kept).

    grid=True uses the energy-grid integral when s.integration says so
    (meant for the mean / most probable curve: it is slower).
    sigma: optional per-point uncertainty (used if s.weight_by_sigma).
    """
    V, I = np.asarray(V, float), np.asarray(I, float)
    ok = np.isfinite(V) & np.isfinite(I)
    if sigma is not None:
        sigma = np.asarray(sigma, float)
        ok &= np.isfinite(sigma) & (sigma > 0)
    if s.fit_v_max is not None:
        ok &= np.abs(V) <= s.fit_v_max
    V, I = V[ok], I[ok]
    if len(V) < 8 or np.ptp(V) <= 0:
        return SLMFit()
    scale = np.max(np.abs(I)) or 1.0
    wts = np.ones(len(V))
    if sigma is not None and s.weight_by_sigma:
        sg = sigma[ok]
        wts = np.median(sg) / sg
    n0 = s.n_molecules
    g = low_bias_conductance(V, I) / n0
    lm = s.fit_method == "Levenberg-Marquardt"
    names = ["eps0", "Gamma"] + (["a"] if s.fit_asymmetry else []) + (["n"] if s.fit_n else [])

    def unpack(p):
        d = dict(zip(names, p))
        return d["eps0"], d["Gamma"], d.get("a", 0.0), d.get("n", n0)

    def resid(p):
        e, G_, a, n = unpack(p)
        return wts * (model_current(V, e, G_, a, n, s, grid) - I) / scale

    lo = {"eps0": 0.0, "Gamma": 1e-6, "a": -0.5, "n": 1e-6}
    hi = {"eps0": 10.0, "Gamma": 5.0, "a": 0.5, "n": 1e16}
    vmax = np.max(np.abs(V))
    starts = []
    if grid and s.integration == "energy grid":
        # the closed form gives the same current (to ~1e-9), so it finds the start and
        # the slow grid integral only polishes it
        quick = fit_slm(V, I, s, sigma=None if sigma is None else sigma[ok])
        if quick.ok:
            starts = [{"eps0": quick.eps0, "Gamma": quick.Gamma, "a": quick.a, "n": quick.n}]
    best = None
    for e0 in ((0.4 * vmax, 0.8 * vmax, 1.5 * vmax, 3.0 * vmax) if not starts else (None,)):
        if e0 is None:
            pass
        elif s.fit_n:
            # start from Gamma = e0 / 10 and let N carry the magnitude
            G_0 = e0 / 10
            n_0 = max(abs(np.polyfit(V, I, 1)[0]) / G0 / slm_conductance(e0, G_0), 1e-6) \
                if np.ptp(V) > 0 else n0
        else:
            gg = g if np.isfinite(g) and 0 < g < 0.9 else 1e-3
            G_0, n_0 = min(max(e0 * np.sqrt(gg / (1 - gg)), 1.1e-6), 4.9), n0
        start = {"eps0": min(e0, 9.9), "Gamma": G_0, "a": 0.0, "n": n_0} if e0 is not None else starts[0]
        p0 = [min(max(start[k], lo[k] * 1.0001 + 1e-12), hi[k] * 0.9999) for k in names]
        try:
            if lm:
                r = least_squares(resid, p0, method="lm", x_scale="jac")
            else:
                r = least_squares(resid, p0, bounds=([lo[k] for k in names], [hi[k] for k in names]),
                                  x_scale="jac")
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
        err = dict(zip(names, np.sqrt(np.clip(np.diag(cov), 0, None))))
    except np.linalg.LinAlgError:
        err = {k: np.nan for k in names}
    ss_tot = np.sum((wts * (I - I.mean()) / scale) ** 2)
    r2 = 1 - 2 * best.cost / ss_tot if ss_tot > 0 else np.nan
    e, G_, a, n = unpack(p)
    return SLMFit(abs(float(e)), abs(float(G_)), float(a), float(n), float(err["eps0"]), float(err["Gamma"]),
                  float(err.get("a", np.nan)), float(r2), True, float(err.get("n", np.nan)))


def most_probable_curve(curves, s: IVModelSettings):
    """Most probable I in each bias bin (Gaussian fit of the I distribution), as in
    Rashid et al. 2025 (SI Fig. S22): returns bias, most probable I, sigma, counts."""
    if not curves:
        return None
    V = np.concatenate([c.V for c in curves])
    I = np.concatenate([c.I for c in curves])
    ok = np.isfinite(V) & np.isfinite(I)
    V, I = V[ok], I[ok]
    edges = np.linspace(V.min(), V.max(), int(s.mp_v_bins) + 1)
    k = np.clip(np.digitize(V, edges) - 1, 0, len(edges) - 2)
    vc, mu, sd, cnt = [], [], [], []
    for j in range(len(edges) - 1):
        ij = I[k == j]
        if len(ij) < 10:
            continue
        lo, hi = np.percentile(ij, [1, 99])
        if hi <= lo:
            continue
        h, e = np.histogram(ij, np.linspace(lo, hi, int(s.mp_i_bins) + 1))
        c = 0.5 * (e[1:] + e[:-1])
        m, sg = c[np.argmax(h)], np.std(ij)
        try:
            from scipy.optimize import curve_fit
            popt, _ = curve_fit(lambda x, A, m_, s_: A * np.exp(-(x - m_) ** 2 / (2 * s_ ** 2)), c, h,
                                p0=[h.max(), m, max(sg, 1e-30)], maxfev=2000)
            if lo <= popt[1] <= hi:
                m, sg = popt[1], abs(popt[2])
        except (RuntimeError, ValueError):
            pass
        vc.append(0.5 * (edges[j] + edges[j + 1]))
        mu.append(m)
        sd.append(sg)
        cnt.append(len(ij))
    if len(vc) < 8:
        return None
    return np.array(vc), np.array(mu), np.array(sd), np.array(cnt)


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
    mean: dict                 # "forward" / "backward" / "most probable" -> (V, I, SLMFit, v_plus, v_minus)
    settings: IVModelSettings
    most_probable: tuple | None = None   # (V, I, sigma, counts) of the most probable curve

    def column(self, name: str) -> np.ndarray:
        return np.array([getattr(f, name) for f in self.fits], float)

    def summary_lines(self) -> list[str]:
        lines = []
        for d, (_, _, f, vp, vm) in self.mean.items():
            e, a = eps0_from_vt(vp, vm)
            lines.append(f"{'' if d == 'most probable' else 'mean '}{d} curve: SLM eps0 = {f.eps0:.3f} +/- {f.eps0_err:.3f} eV, "
                         f"Gamma = {1e3 * f.Gamma:.3g} +/- {1e3 * f.Gamma_err:.2g} meV, a = {f.a:.3f}, "
                         + (f"N = {f.n:.3g}, " if self.settings.fit_n else "") +
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
        fits.append(SLMFit())
    single = IVModelSettings(**{**asdict(s), "integration": "analytic"})
    for k, c in enumerate(curves[:s.max_curves]):
        fits[k] = fit_slm(c.V, c.I, single)
    mean = {}
    if mean_fn is not None:
        for d, (vg, m, _) in mean_fn(curves).items():
            p, q = transition_voltages(vg, m, s.tvs_v_min)
            mean[d] = (vg, m, fit_slm(vg, m, s, grid=True), p, q)
    mp = most_probable_curve(curves, s)
    if mp is not None:
        vg, m, sd, _ = mp
        p, q = transition_voltages(vg, m, s.tvs_v_min)
        mean["most probable"] = (vg, m, fit_slm(vg, m, s, sigma=sd, grid=True), p, q)
        most_probable = mp
    else:
        most_probable = None
    return IVModelResult(np.array([c.index for c in curves]), np.array([c.direction for c in curves]),
                         vp, vm, et, at, gl, fits, mean, s, most_probable)
