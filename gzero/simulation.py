"""Fake break-junction data with known answers, for testing.

One cycle: approach until the contact snaps in, pull through the atomic
steps (sometimes with a molecular plateau), hold with 1/f noise where
NP ~ G^n, retract, short hold. The result has a piezo channel, so it goes
through the normal trace detection.
"""

from __future__ import annotations

from dataclasses import dataclass, asdict

import numpy as np

from . import G0
from .conversion import Recording


@dataclass
class SimulationSettings:
    n_cycles: int = 200
    fs: float = 10000.0
    speed: float = 2.0                 # junction nm / s
    displacement_ratio: float = 1.0    # junction nm per piezo nm
    beta: float = 10.0                 # tunnelling decay (1/nm, ln units)
    contact_atoms: tuple = (1, 4)      # random contact size range (G0)
    atom_step: float = 0.08            # nm per atomic plateau
    snap_back: float = 0.3             # nm gap opened at rupture
    p_molecule: float = 0.7            # fraction of traces with a molecular plateau
    molecule_logG: float = -3.5        # centre of molecular conductance (log G/G0)
    molecule_sigma: float = 0.25       # spread (decades)
    plateau_length: float = 0.4        # mean molecular plateau length (nm)
    plateau_length_sd: float = 0.15
    plateau_slope: float = 0.3         # decades / nm drop along the plateau
    second_molecule_logG: float | None = None   # optional second class (clustering test)
    p_second: float = 0.0
    hold_s: float = 1.0                # hold after the pull (0 = no hold)
    pull_distance: float = 1.0         # nm pulled before the hold
    retract_distance: float = 1.0      # extra nm after the hold
    noise_exponent: float = 1.5        # NP ~ G^n during holds
    noise_amplitude: float = 0.05      # relative rms fluctuation at G = 10^molecule_logG
    noise_floor_logG: float = -5.5
    floor_noise: float = 0.3           # relative noise of the floor
    measurement_noise: float = 0.02    # relative white noise on G
    bias: float = 0.1
    seed: int = 1

    def to_dict(self) -> dict:
        return asdict(self)


def _pink(n: int, rng) -> np.ndarray:
    """Unit-variance 1/f noise."""
    if n < 4:
        return rng.standard_normal(n)
    f = np.fft.rfftfreq(n)
    spec = rng.standard_normal(len(f)) + 1j * rng.standard_normal(len(f))
    f[0] = f[1]
    spec /= np.sqrt(f)
    x = np.fft.irfft(spec, n)
    x -= x.mean()
    return x / (x.std() or 1.0)


def simulate(s: SimulationSettings) -> Recording:
    rng = np.random.default_rng(s.seed)
    dt = 1.0 / s.fs
    step = s.speed * dt                     # junction nm per sample
    floor = 10.0 ** s.noise_floor_logG
    G_parts, z_parts = [], []
    z_now = 0.0                             # junction coordinate (nm), larger = more open
    truth = {"molecule_logG": [], "plateau_length": [], "hold_logG": [], "cycle_kind": []}

    def tunnel(gap):
        return np.exp(-s.beta * np.clip(gap, 0, None))

    for _ in range(s.n_cycles):
        # 1. approach: gap closes from current opening to contact
        gap0 = s.pull_distance + s.retract_distance
        n_ap = int(gap0 / step)
        gap = gap0 - step * np.arange(n_ap)
        Ga = tunnel(gap)
        snap = np.flatnonzero(Ga > 0.1)
        n_atoms = int(rng.integers(s.contact_atoms[0], s.contact_atoms[1] + 1))
        if len(snap):
            Ga[snap[0]:] = n_atoms
        G_parts.append(Ga)
        z_parts.append(z_now - step * np.arange(n_ap))
        z_now -= step * n_ap

        # 2. pull: atomic steps, rupture, optional molecule, tunnelling
        n_pull = int(s.pull_distance / step)
        d = step * np.arange(n_pull)
        Gp = np.empty(n_pull)
        x_rupt = n_atoms * s.atom_step * rng.uniform(0.7, 1.3)
        steps_edges = np.linspace(0, x_rupt, n_atoms + 1)
        for k in range(n_atoms):
            sel = (d >= steps_edges[k]) & (d < steps_edges[k + 1])
            Gp[sel] = (n_atoms - k) * rng.uniform(0.95, 1.02)
        after = d >= x_rupt
        r = rng.random()
        if r < s.p_molecule:
            if s.second_molecule_logG is not None and rng.random() < s.p_second:
                g_mol = s.second_molecule_logG + s.molecule_sigma * rng.standard_normal()
                kind = "molecule2"
            else:
                g_mol = s.molecule_logG + s.molecule_sigma * rng.standard_normal()
                kind = "molecule"
            L = abs(s.plateau_length + s.plateau_length_sd * rng.standard_normal())
            on_plat = after & (d < x_rupt + L)
            Gp[on_plat] = 10.0 ** (g_mol - s.plateau_slope * (d[on_plat] - x_rupt))
            past = after & (d >= x_rupt + L)
            g_end = 10.0 ** (g_mol - s.plateau_slope * L)
            Gp[past] = g_end * tunnel(s.snap_back * 0.3 + d[past] - x_rupt - L)
            truth["molecule_logG"].append(g_mol)
            truth["plateau_length"].append(L)
        else:
            kind = "tunnel"
            Gp[after] = tunnel(s.snap_back + d[after] - x_rupt)
        truth["cycle_kind"].append(kind)
        G_parts.append(Gp)
        z_parts.append(z_now + d)
        z_now += step * n_pull

        # 3. hold: 1/f noise with NP ~ G^n
        if s.hold_s > 0:
            n_h = int(s.hold_s * s.fs)
            g_hold = max(Gp[-1], floor)
            lg = np.log10(g_hold)
            ref = 10.0 ** s.molecule_logG
            # rms fluctuation so that variance = (amp*ref)^2 * (G/ref)^n
            rms = s.noise_amplitude * ref * (g_hold / ref) ** (s.noise_exponent / 2.0)
            rms *= np.exp(0.25 * rng.standard_normal())    # junction-to-junction scatter
            Gh = g_hold + rms * _pink(n_h, rng)
            G_parts.append(np.clip(Gh, g_hold * 0.05, None))
            z_parts.append(np.full(n_h, z_now))
            truth["hold_logG"].append(lg)

        # 4. retract
        n_re = int(s.retract_distance / step)
        d2 = step * np.arange(1, n_re + 1)
        G_parts.append(max(Gp[-1], floor) * tunnel(d2))
        z_parts.append(z_now + d2)
        z_now += step * n_re

        # 5. short hold
        n_h2 = int(0.05 * s.fs)
        G_parts.append(np.zeros(n_h2))
        z_parts.append(np.full(n_h2, z_now))

    G = np.concatenate(G_parts)
    z = np.concatenate(z_parts)
    G = G * (1.0 + s.measurement_noise * rng.standard_normal(len(G)))
    G = G + floor * (1.0 + s.floor_noise * rng.standard_normal(len(G)))
    G = np.abs(G)
    piezo = z / s.displacement_ratio
    rec = Recording(name="simulation", fs=s.fs, G=G, piezo=piezo, bias=s.bias,
                    current=G * G0 * s.bias,
                    meta={"simulation": s.to_dict(),
                          "truth": {k: np.asarray(v) if k != "cycle_kind" else v for k, v in truth.items()}})
    return rec
