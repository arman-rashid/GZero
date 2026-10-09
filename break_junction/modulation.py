"""Piezo modulation: lock-in of ln G and z at the modulation frequency.

beta = -d ln G / dz. For tunnelling that's about 2 kappa, on a molecular
plateau it's small.
"""

from __future__ import annotations

from dataclasses import dataclass, asdict

import numpy as np

from .conversion import Recording


@dataclass
class ModulationSettings:
    f_mod: float = 0.0          # Hz; 0 = detect from the piezo spectrum
    periods_per_window: int = 4
    f_search_min: float = 10.0  # Hz, search range for automatic detection
    g_min: float = -6.0         # only report windows with mean log G in range
    g_max: float = 1.0
    displacement_ratio: float = 1.0   # junction nm per piezo nm

    def to_dict(self) -> dict:
        return asdict(self)


def spectrum(x: np.ndarray, fs: float):
    x = np.asarray(x, float) - np.mean(x)
    w = np.hanning(len(x))
    F = np.fft.rfft(x * w)
    f = np.fft.rfftfreq(len(x), 1.0 / fs)
    amp = 2.0 * np.abs(F) / np.sum(w)
    return f, amp


def detect_frequency(piezo: np.ndarray, fs: float, f_min: float = 10.0,
                     min_peak_ratio: float = 20.0) -> float:
    """Strongest spectral line above f_min; 0 if there is no clear line.

    The piezo ramps and holds give a dense comb of harmonics that falls off
    with frequency, so a real modulation must stand out min_peak_ratio
    times above the spectrum in its own neighbourhood (not just globally).
    """
    n = min(len(piezo), int(fs * 2))  # 2 s is plenty
    f, a = spectrum(np.diff(piezo[:n]), fs)   # differentiate: suppress the ramps
    idx = np.flatnonzero(f >= f_min)
    if len(idx) < 10:
        return 0.0
    i = idx[np.argmax(a[idx])]
    lo, hi = max(idx[0], i - 60), min(len(a), i + 61)
    neighbours = np.r_[a[lo:max(lo, i - 3)], a[min(hi, i + 4):hi]]
    if len(neighbours) == 0 or a[i] < min_peak_ratio * np.median(neighbours):
        return 0.0
    return float(f[i])


def lockin(rec: Recording, ms: ModulationSettings, start: int = 0, stop: int | None = None):
    """Return dict with per-window time, mean log G, amplitude of z and ln G, beta."""
    if rec.piezo is None:
        raise ValueError("Piezo modulation analysis needs a piezo channel.")
    stop = rec.n if stop is None else stop
    z = np.asarray(rec.piezo[start:stop], float) * ms.displacement_ratio
    lnG = np.log(np.clip(rec.G[start:stop], 1e-12, None))
    fm = ms.f_mod or detect_frequency(z, rec.fs, ms.f_search_min)
    if fm <= 0:
        raise ValueError("No piezo modulation was found in this recording (no clear spectral line). "
                         "Enter the modulation frequency if you know it.")
    per = rec.fs / fm
    w = max(8, int(round(ms.periods_per_window * per)))
    t = np.arange(len(z)) / rec.fs
    c, s = np.cos(2 * np.pi * fm * t), np.sin(2 * np.pi * fm * t)
    out = {"t": [], "logG": [], "amp_z": [], "amp_lnG": [], "beta": [], "phase": []}
    for k in range(len(z) // w):
        sl = slice(k * w, (k + 1) * w)
        zz = z[sl] - np.polyval(np.polyfit(t[sl], z[sl], 1), t[sl])
        gg = lnG[sl] - np.polyval(np.polyfit(t[sl], lnG[sl], 1), t[sl])
        zc = 2 * np.mean(zz * c[sl]) + 2j * np.mean(zz * s[sl])
        gc = 2 * np.mean(gg * c[sl]) + 2j * np.mean(gg * s[sl])
        lg = np.mean(lnG[sl]) / np.log(10)
        if not (ms.g_min <= lg <= ms.g_max) or abs(zc) == 0:
            continue
        ratio = gc / zc
        out["t"].append((start + k * w + w / 2) / rec.fs)
        out["logG"].append(lg)
        out["amp_z"].append(abs(zc))
        out["amp_lnG"].append(abs(gc))
        out["beta"].append(-ratio.real)          # sign: G drops when z increases -> beta > 0
        out["phase"].append(np.degrees(np.angle(ratio)))
    res = {k: np.asarray(v) for k, v in out.items()}
    res["f_mod"] = fm
    return res
