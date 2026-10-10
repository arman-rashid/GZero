"""Per-trace plateau statistics: junction yield and conductance over the run.

A trace "forms a junction" if it stays flat inside the molecular window
for at least min_plateau nm. Flat means the local slope |d log G / dz|
is below max_slope, which separates molecular plateaus (well under
1 decade/nm) from tunnelling through the window (beta/ln10, about
4-9 decades/nm). Without the slope test, a tunnelling-only trace counts
as a short plateau, because it takes about 0.1-0.2 nm to fall through a
typical window.

Per trace this gives the plateau length, the plateau conductance (the
median log G of the flat points) and a yes/no. The yield is the fraction
of traces with a junction (Kamenetska et al., PRL 102, 126803 (2009)).
In blocks of consecutive traces it shows whether the yield and the
conductance drift during the measurement.
"""

from __future__ import annotations

from dataclasses import dataclass, asdict

import numpy as np

from .traces import Trace


@dataclass
class JunctionSettings:
    g_lo: float = -4.5              # molecular window, log G
    g_hi: float = -2.5
    max_slope: float = 2.0          # decades / nm; flatter counts as plateau
    slope_window: float = 0.05      # nm over which the slope is taken
    min_plateau: float = 0.1        # nm of flat plateau needed to count as a junction
    block_size: int = 100           # traces per block for the evolution plots
    g_bins: int = 60                # bins of the plateau conductance histogram

    def to_dict(self) -> dict:
        return asdict(self)


def local_slope(z: np.ndarray, logG: np.ndarray, width: float) -> np.ndarray:
    """|d log G / dz| at every point, from the values width/2 before and after (decades/nm)."""
    z = np.maximum.accumulate(np.asarray(z, float))   # distance must not go back
    if len(z) < 3 or z[-1] - z[0] <= 0:
        return np.full(len(z), np.inf)
    zu, first = np.unique(z, return_index=True)
    gu = np.asarray(logG, float)[first]
    if len(zu) < 3:
        return np.full(len(z), np.inf)
    h = width / 2
    lo = np.interp(z - h, zu, gu)
    hi = np.interp(z + h, zu, gu)
    span = np.minimum(z + h, zu[-1]) - np.maximum(z - h, zu[0])
    with np.errstate(divide="ignore", invalid="ignore"):
        s = np.abs(hi - lo) / span
    return np.where(span > 0.5 * width, s, np.inf)


@dataclass
class JunctionResult:
    index: np.ndarray           # trace numbers
    length: np.ndarray          # flat plateau length in the window (nm)
    logG: np.ndarray            # plateau conductance (median log G of the flat points), NaN = none
    junction: np.ndarray        # bool
    settings: JunctionSettings

    @property
    def n(self) -> int:
        return len(self.index)

    def yield_(self) -> tuple[float, float]:
        """Fraction of traces with a junction and its binomial standard error."""
        if not self.n:
            return np.nan, np.nan
        p = float(np.mean(self.junction))
        return p, float(np.sqrt(p * (1 - p) / self.n))

    def blocks(self):
        """Per block of consecutive traces: centre trace number, yield, its SE, median plateau log G."""
        b = max(1, int(self.settings.block_size))
        out = []
        for s in range(0, self.n, b):
            sl = slice(s, min(self.n, s + b))
            j = self.junction[sl]
            p = float(np.mean(j))
            g = self.logG[sl][j]
            out.append((float(np.mean(self.index[sl])), p, float(np.sqrt(p * (1 - p) / len(j))),
                        float(np.median(g)) if len(g) else np.nan, int(len(j))))
        return np.array(out).reshape(-1, 5)

    def summary_lines(self) -> list[str]:
        p, se = self.yield_()
        g = self.logG[self.junction]
        L = self.length[self.junction]
        s = self.settings
        lines = [f"{self.n} traces, {int(self.junction.sum())} with a plateau of at least {s.min_plateau:g} nm "
                 f"in [{s.g_lo:g}, {s.g_hi:g}] (slope below {s.max_slope:g} dec/nm)",
                 f"junction yield = {100 * p:.1f} +/- {100 * se:.1f} %"]
        if len(g):
            lines.append(f"plateau log G: median {np.median(g):.3f}, mean {np.mean(g):.3f}, "
                         f"SD {np.std(g):.3f} (G = {10 ** np.median(g):.3g} G0)")
            lines.append(f"plateau length: median {np.median(L):.3f} nm, mean {np.mean(L):.3f} nm")
        bl = self.blocks()
        if len(bl) >= 3:
            ok = np.isfinite(bl[:, 3])
            if ok.sum() >= 3:
                drift = np.polyfit(bl[ok, 0], bl[ok, 3], 1)[0] * 1000
                lines.append(f"conductance drift: {drift:+.3f} decades per 1000 traces")
            lines.append(f"yield in blocks of {s.block_size}: {100 * bl[:, 1].min():.0f} - "
                         f"{100 * bl[:, 1].max():.0f} %")
        return lines


def analyse(traces: list[Trace], s: JunctionSettings) -> JunctionResult:
    n = len(traces)
    length = np.zeros(n)
    logg = np.full(n, np.nan)
    for i, t in enumerate(traces):
        y, z = t.logG, t.z
        if len(y) < 5:
            continue
        flat = (y >= s.g_lo) & (y <= s.g_hi) & (local_slope(z, y, s.slope_window) <= s.max_slope)
        if not flat.any():
            continue
        dz = np.abs(np.diff(z))
        step = np.r_[dz, dz[-1:]] if len(dz) else np.zeros(1)
        length[i] = float(np.sum(step[flat]))
        logg[i] = float(np.median(y[flat]))
    junction = length >= s.min_plateau
    logg[~junction] = np.nan
    return JunctionResult(np.array([t.index for t in traces]), length, logg, junction, s)
