"""I-V sweeps.

The bias is split at its turning points; each monotonic piece is a forward
or backward curve. dI/dV comes from a Savitzky-Golay derivative.
"""

from __future__ import annotations

from dataclasses import dataclass, asdict, field

import numpy as np
from scipy.signal import find_peaks, savgol_filter

from . import G0
from .conversion import Recording


@dataclass
class IVSettings:
    min_amplitude: float = 0.05        # V; turning points must be this far apart
    smooth_points: int = 11            # Savitzky-Golay window for current / dI/dV
    poly_order: int = 2
    v_min: float | None = None         # keep curves covering at least this range
    v_max: float | None = None
    v_bins: int = 100
    i_bins: int = 100
    log_current: bool = True

    def to_dict(self) -> dict:
        return asdict(self)


@dataclass
class IVCurve:
    index: int
    direction: str          # "forward" | "backward"
    V: np.ndarray
    I: np.ndarray
    dIdV: np.ndarray = field(default=None)

    @property
    def G(self) -> np.ndarray:
        with np.errstate(divide="ignore", invalid="ignore"):
            return self.I / self.V / G0


def split_sweeps(V: np.ndarray, I: np.ndarray, s: IVSettings) -> list[IVCurve]:
    V = np.asarray(V, float)
    I = np.asarray(I, float)
    if np.ptp(V) < s.min_amplitude:
        return []
    width = max(3, s.smooth_points)
    Vs = savgol_filter(V, width | 1, min(s.poly_order, (width | 1) - 1)) if len(V) > width else V
    pk, _ = find_peaks(Vs, prominence=s.min_amplitude)
    vl, _ = find_peaks(-Vs, prominence=s.min_amplitude)
    turns = np.unique(np.r_[0, pk, vl, len(V) - 1])
    curves = []
    for a, b in zip(turns[:-1], turns[1:]):
        if b - a < max(10, width):
            continue
        v, i = V[a:b + 1], I[a:b + 1]
        direction = "forward" if v[-1] > v[0] else "backward"
        if s.v_min is not None and v.min() > s.v_min:
            continue
        if s.v_max is not None and v.max() < s.v_max:
            continue
        w = min(len(v) - (1 - len(v) % 2), width | 1)
        if w >= 5:
            i_s = savgol_filter(i, w, min(s.poly_order, w - 1))
            di = savgol_filter(i, w, min(s.poly_order, w - 1), deriv=1)
            dv = savgol_filter(v, w, min(s.poly_order, w - 1), deriv=1)
            with np.errstate(divide="ignore", invalid="ignore"):
                didv = np.where(np.abs(dv) > 0, di / dv, np.nan)
        else:
            i_s, didv = i, np.gradient(i, v)
        curves.append(IVCurve(len(curves), direction, v, i_s, didv))
    return curves


def curves_from_recording(rec: Recording, s: IVSettings) -> list[IVCurve]:
    if rec.current is None or np.ndim(rec.bias) == 0:
        raise ValueError("I-V analysis needs a time-resolved bias and current "
                         "(a TDMS file, or a text file with a bias column).")
    return split_sweeps(rec.bias_array(), rec.current, s)


def iv_hist2d(curves: list[IVCurve], s: IVSettings, quantity: str = "I"):
    """2D histogram of I (or log|I|, G, dI/dV) versus V over all curves."""
    if not curves:
        return None
    V = np.concatenate([c.V for c in curves])
    if quantity == "G":
        Y = np.log10(np.clip(np.abs(np.concatenate([c.G for c in curves])), 1e-12, None))
    elif quantity == "dI/dV":
        Y = np.concatenate([c.dIdV for c in curves]) / G0
        Y = np.log10(np.clip(np.abs(Y), 1e-12, None))
    else:
        Y = np.concatenate([c.I for c in curves])
        if s.log_current:
            Y = np.log10(np.clip(np.abs(Y), 1e-15, None))
    ok = np.isfinite(V) & np.isfinite(Y)
    H, ve, ye = np.histogram2d(V[ok], Y[ok], [s.v_bins, s.i_bins])
    return H / len(curves), ve, ye


def mean_curve(curves: list[IVCurve], n: int = 200):
    """Average I(V) on a common voltage grid (forward and backward separately)."""
    out = {}
    for d in ("forward", "backward"):
        cs = [c for c in curves if c.direction == d]
        if not cs:
            continue
        lo = max(c.V.min() for c in cs)
        hi = min(c.V.max() for c in cs)
        if hi <= lo:
            continue
        vg = np.linspace(lo, hi, n)
        I = np.array([np.interp(vg, np.sort(c.V), c.I[np.argsort(c.V)]) for c in cs])
        out[d] = (vg, I.mean(axis=0), I.std(axis=0))
    return out
