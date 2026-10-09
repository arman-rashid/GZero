"""Raw data -> conductance.

Calibrate the current channels (linear gain or the diode polynomials),
merge the low- and high-gain stages of a dual-stage amplifier (or take the
single-stage current), correct for the series resistor
(G = I / (E - I R) / G0), filter if wanted. Quantities a file already
contains (conductance, bias voltage, piezo, current, time) are used as they
are; only the missing ones are calculated. Recording is the format the
rest of the code works with.
"""

from __future__ import annotations

import os
import re
from dataclasses import dataclass, field, asdict

import numpy as np
from scipy.ndimage import median_filter, uniform_filter1d
from scipy.signal import savgol_filter

from . import G0
from .bj_io import TdmsData


# Common in-memory representation used by every analysis module

@dataclass
class Recording:
    """One continuous measurement converted to conductance.

    G is in units of G0 (linear). piezo is the piezo position in nm
    (None if the file has no piezo column; then distance is derived from
    the sample index and the "piezo rate").
    """

    name: str
    fs: float
    G: np.ndarray
    piezo: np.ndarray | None = None
    bias: np.ndarray | float = 0.1          # applied bias (V)
    current: np.ndarray | None = None        # junction current (A)
    meta: dict = field(default_factory=dict)
    events: list = field(default_factory=list)

    @property
    def n(self) -> int:
        return len(self.G)

    @property
    def time(self) -> np.ndarray:
        return np.arange(self.n) / self.fs

    def logG(self, floor: float = 1e-10) -> np.ndarray:
        return np.log10(np.clip(self.G, floor, None))

    def bias_array(self) -> np.ndarray:
        if np.ndim(self.bias) == 0:
            return np.full(self.n, float(self.bias))
        return np.asarray(self.bias, float)


# Filters (Savitzky-Golay / median, optionally applied to log|x|)

@dataclass
class FilterSettings:
    kind: str = "none"        # "none" | "savgol" | "median" | "moving average"
    side_points: int = 5      # window = 2*side_points + 1
    poly_order: int = 2
    in_log: bool = False      # filter log10|x| then transform back (keeps sign)


def apply_filter(x: np.ndarray, fs: FilterSettings) -> np.ndarray:
    if fs is None or fs.kind == "none" or fs.side_points < 1:
        return x
    x = np.asarray(x, float)
    win = 2 * int(fs.side_points) + 1
    if win > len(x):
        return x
    sign = None
    if fs.in_log:
        sign = np.sign(x)
        sign[sign == 0] = 1
        y = np.log10(np.clip(np.abs(x), 1e-30, None))
    else:
        y = x
    if fs.kind == "savgol":
        order = min(int(fs.poly_order), win - 1)
        y = savgol_filter(y, win, order, mode="interp")
    elif fs.kind == "median":
        y = median_filter(y, size=win, mode="nearest")
    elif fs.kind == "moving average":
        y = uniform_filter1d(y, win, mode="nearest")
    else:
        raise ValueError(f"unknown filter {fs.kind!r}")
    if fs.in_log:
        y = sign * 10.0 ** y
    return y


# Diode / logarithmic amplifier calibration (module 2 of the EC4 setup)

@dataclass
class DiodeCalibration:
    """Piece-wise polynomial V -> I calibration stored in Calibration_M2.

    Regions on u = V - V_offset:  |u| < t_low -> middle polynomial gives
    I directly (A); t_low <= |u| < t_high -> low polynomial gives log10|I|;
    |u| >= t_high -> high polynomial gives log10|I|. Negative u uses the
    negative-side polynomials evaluated at |u| and returns -I.
    Coefficients are in ascending order (c0 + c1 u + ...), as stored.
    """

    pos_high: list
    pos_low: list
    middle: list
    neg_low: list
    neg_high: list
    v_offset: float = 0.0
    t_low: float = 0.07
    t_high: float = 0.5
    invert: bool = False

    @staticmethod
    def from_tdms(cal: dict, prefix: str) -> "DiodeCalibration | None":
        def get(key):
            v = cal.get(f"{prefix}-{key}", [])
            return v if isinstance(v, list) else []

        parts = [get(k) for k in ("PH", "PL", "M", "NL", "NH")]
        if not all(parts):
            return None
        thr = get("Threshold") or [0.07, 0.5]
        off = get("V_offset") or [0.0]
        return DiodeCalibration(*parts, v_offset=off[0], t_low=thr[0], t_high=thr[1])

    def __call__(self, v: np.ndarray) -> np.ndarray:
        u = np.asarray(v, float) - self.v_offset
        a = np.abs(u)
        out = np.empty_like(u)
        poly = np.polynomial.polynomial.polyval

        mid = a < self.t_low
        out[mid] = poly(u[mid], self.middle)
        for sgn, sel_sign, low, high in ((1.0, u >= 0, self.pos_low, self.pos_high),
                                         (-1.0, u < 0, self.neg_low, self.neg_high)):
            lo = sel_sign & ~mid & (a < self.t_high)
            hi = sel_sign & (a >= self.t_high)
            out[lo] = sgn * 10.0 ** poly(a[lo], low)
            out[hi] = sgn * 10.0 ** poly(a[hi], high)
        return -out if self.invert else out


# Amplifier currents -> conductance (shared by TDMS and column files)

AMPLIFIERS = ["dual stage", "single stage"]
CONDUCTANCE_UNITS = ["G/G0", "S", "log G/G0"]


def currents_to_conductance(ia, ib, E, R=0.0, amplifier="dual stage", combine="auto", b_ok=None,
                            switch_logG=None, a_saturated=None, g_ceiling=100.0):
    """Conductance (G0) from the amplifier current(s) in A.

    ia is the low-gain stage Ia (or the only current of a single-stage
    amplifier), ib the high-gain stage Ib. Dual stage with combine "auto":
    Ib wherever b_ok (Ib not saturated) and, if switch_logG is set, only
    below that log G; Ia elsewhere. Single stage uses Ia alone. E is the
    bias (scalar or array), R the series resistance. Samples where Ia
    saturates (a_saturated) or where R carries the whole bias get
    g_ceiling. Returns (G, I, fraction of samples taken from Ib or None).
    """
    def conductance(I):
        with np.errstate(divide="ignore", invalid="ignore"):
            return I / (E - I * R) / G0

    if ia is None and ib is None:
        raise ValueError("No current channel.")
    frac_b = None
    if amplifier == "single stage" or ib is None or combine == "a":
        I = ia if ia is not None else ib
        uses_a = ia is not None
    elif ia is None or combine == "b":
        I, uses_a = ib, False
    else:
        use_b = np.ones(len(ib), bool) if b_ok is None else np.asarray(b_ok, bool).copy()
        if switch_logG is not None:
            ga = np.abs(conductance(ia))
            use_b &= np.log10(np.clip(ga, 1e-12, None)) < switch_logG
        I = np.where(use_b, ib, ia)
        frac_b = float(np.mean(use_b))
        uses_a = True

    G = conductance(I)
    # Where E - I*R <= 0 the limiting resistor carries the whole bias (closed
    # contact beyond what can be resolved): clamp to the ceiling.
    saturated = ~np.isfinite(G) | ((E - I * R <= 0) & (I * E > 0))
    if a_saturated is not None and uses_a:
        # low-gain stage saturated: the contact conducts more than it can measure
        saturated |= np.asarray(a_saturated, bool)
    G = np.where(saturated, g_ceiling, np.abs(G))
    return np.clip(G, None, g_ceiling), I, frac_b


def conductance_from_unit(x: np.ndarray, unit: str) -> np.ndarray:
    """Conductance column / channel in the given unit -> G/G0."""
    unit = unit.strip()
    if unit in ("G/G0", "G0"):
        return np.abs(x)
    if unit in ("log G/G0", "logG"):
        return 10.0 ** x
    if unit in ("S", "siemens", ""):
        return np.abs(x) / G0
    raise ValueError(f"unknown conductance unit {unit!r}")


# TDMS -> Recording

@dataclass
class ConversionSettings:
    module: int = 1                    # which measurement module (M_E{m}, M_i{m}a/b)
    amplifier: str = "dual stage"      # "dual stage" (Ia low + Ib high gain) | "single stage" (Ia only)
    # channel names; "auto" = the EC4 names of the module, "none" = not used
    current_a_channel: str = "auto"    # M_i{m}a
    current_b_channel: str = "auto"    # M_i{m}b (dual stage only)
    conductance_channel: str = "auto"  # auto: "G" only if the file has no current channels
    conductance_unit: str = "auto"     # "auto" (channel unit, else S) | "S" | "G/G0" | "log G/G0"
    piezo_channel: str = "auto"        # U_Piezo
    calibration: str = "linear"        # "linear" | "diode"
    # linear calibration: I = (V - offset) * gain * corr
    gain_a: float | None = None        # A/V; None -> channel scale_factor from file
    gain_b: float | None = None
    offset_a: float = 0.0
    offset_b: float = 0.0
    corr_a: float = 1.0
    corr_b: float = 1.0
    amplification: float = 1.0         # extra overall factor
    invert_current: bool = False
    # channel combination (dual stage)
    combine: str = "auto"              # "auto" | "a" | "b"
    b_saturation_V: float = 10.0       # |V_b| above this -> high-gain channel saturated
    a_saturation_V: float | None = 10.4  # |V_a| above this -> contact beyond range -> g_ceiling
    switch_logG: float | None = None   # optional: use b only below this log10(G/G0)
    # bias
    bias_source: str = "measured"      # "measured" | "fixed"
    bias_channel: str = "auto"         # "auto" = M_E1-E2 if present, else M_E{module}; or a channel name
    fixed_bias: float = 0.1
    bias_filter: FilterSettings = field(default_factory=lambda: FilterSettings("moving average", 25))
    # series (current-limiting) resistance compensation
    series_compensation: bool = True
    series_resistance: float | None = None   # ohm; None -> from Setup group
    g_ceiling: float = 100.0                 # G/G0 assigned to saturated samples
    # piezo
    piezo_scale: float | None = None   # nm/V; None -> channel scale_factor
    # filters
    current_filter: FilterSettings = field(default_factory=FilterSettings)
    conductance_filter: FilterSettings = field(default_factory=FilterSettings)

    def to_dict(self) -> dict:
        return asdict(self)

    @staticmethod
    def from_dict(d: dict) -> "ConversionSettings":
        d = dict(d)
        for k in ("bias_filter", "current_filter", "conductance_filter"):
            if isinstance(d.get(k), dict):
                d[k] = FilterSettings(**d[k])
        return ConversionSettings(**d)


def _scale(props: dict, default: float = 1.0) -> float:
    try:
        return float(props.get("scale_factor", default))
    except (TypeError, ValueError):
        return default


def _channel(ch: dict, setting: str, auto: str | None) -> str | None:
    """Channel for a setting: "auto" -> auto (if present), "none"/"" -> None, else the name."""
    s = (setting or "").strip()
    if s.lower() in ("", "none", "off"):
        return None
    if s.lower() == "auto":
        return auto if auto in ch else None
    if s not in ch:
        raise KeyError(f"Channel '{s}' not found; available: {sorted(ch)}")
    return s


def _tdms_piezo(td: TdmsData, cs: ConversionSettings):
    name = _channel(td.channels, cs.piezo_channel, "U_Piezo")
    if name is None:
        return None
    k = cs.piezo_scale if cs.piezo_scale is not None else _scale(td.channel_props[name], 1.0)
    return td.channels[name] * k


def _tdms_bias(td: TdmsData, cs: ConversionSettings, n: int):
    """(bias array, channel name or "fixed"). A measured bias channel replaces the fixed value."""
    ch = td.channels
    if cs.bias_source == "measured":
        if cs.bias_channel.strip().lower() not in ("", "auto"):
            name_e = _channel(ch, cs.bias_channel, None)
        else:
            # the EC4 software itself computes G from the differential bias E1-E2
            name_e = "M_E1-E2" if "M_E1-E2" in ch else f"M_E{cs.module}"
        if name_e in ch:
            return apply_filter(ch[name_e], cs.bias_filter), name_e
    return np.full(n, float(cs.fixed_bias)), "fixed"


def tdms_to_recording(td: TdmsData, cs: ConversionSettings) -> Recording:
    """Convert raw TDMS channels to a :class:`Recording` (G in G0)."""
    m = cs.module
    ch = td.channels
    g_auto = cs.conductance_channel.strip().lower() == "auto"
    name_g = None if g_auto else _channel(ch, cs.conductance_channel, None)
    name_a = _channel(ch, cs.current_a_channel, f"M_i{m}a")
    name_b = _channel(ch, cs.current_b_channel, f"M_i{m}b") if cs.amplifier == "dual stage" else None
    if name_g is None and name_a is None and name_b is None:
        if g_auto and "G" in ch:
            name_g = "G"
        else:
            raise KeyError(f"Module {m} current channels not found; available: {sorted(ch)}")
    if name_g is not None:
        return _precomputed_g_recording(td, cs, name_g, name_a or name_b)

    va = ch[name_a] if name_a else None
    vb = ch[name_b] if name_b else None

    # calibration: raw volts -> amperes
    ia = ib = None
    if cs.calibration == "diode":
        cal = td.calibration.get(m, {})
        if va is not None:
            dc = DiodeCalibration.from_tdms(cal, f"Modula_i{m}")
            ia = dc(va) if dc else None
        if vb is not None:
            for r in range(4):
                dc = DiodeCalibration.from_tdms(cal, f"Modula_i{m}b{r}")
                if dc:
                    ib = dc(vb)
                    break
        if ia is None and ib is None:
            raise ValueError(f"No diode calibration polynomials for module {m} in the file.")
    else:
        if va is not None:
            g = cs.gain_a if cs.gain_a is not None else _scale(td.channel_props[name_a], 1e-7)
            ia = (va - cs.offset_a) * g * cs.corr_a
        if vb is not None:
            g = cs.gain_b if cs.gain_b is not None else _scale(td.channel_props[name_b], 1e-10)
            ib = (vb - cs.offset_b) * g * cs.corr_b

    sign = -1.0 if cs.invert_current else 1.0
    if ia is not None:
        ia = sign * cs.amplification * apply_filter(ia, cs.current_filter)
    if ib is not None:
        ib = sign * cs.amplification * apply_filter(ib, cs.current_filter)

    n = len(va if va is not None else vb)
    E, bias_name = _tdms_bias(td, cs, n)

    R = cs.series_resistance if cs.series_resistance is not None else td.series_resistance()
    if not cs.series_compensation:
        R = 0.0

    linear = cs.calibration != "diode"
    b_ok = np.abs(vb) < cs.b_saturation_V if (vb is not None and linear) else None
    a_sat = None
    if cs.a_saturation_V is not None and va is not None and linear:
        a_sat = np.abs(va) >= cs.a_saturation_V
    G, I, used_b = currents_to_conductance(
        ia, ib, E, R, amplifier=cs.amplifier, combine=cs.combine, b_ok=b_ok,
        switch_logG=cs.switch_logG, a_saturated=a_sat, g_ceiling=cs.g_ceiling)
    G = apply_filter(G, cs.conductance_filter)

    meta = {
        "source": td.path,
        "series_resistance_ohm": R,
        "module": m,
        "amplifier": cs.amplifier,
        "current_channels": [c for c in (name_a, name_b) if c],
        "bias_channel": bias_name,
        "fraction_high_gain": used_b,
        "file_properties": {k: str(v) for k, v in td.file_props.items()},
        "setup": td.setup,
        "conversion": cs.to_dict(),
    }
    return Recording(
        name=os.path.splitext(os.path.basename(td.path))[0],
        fs=td.fs, G=G, piezo=_tdms_piezo(td, cs), bias=E, current=I, meta=meta, events=td.events,
    )


def _precomputed_g_recording(td: TdmsData, cs: ConversionSettings, name_g: str = "G",
                             name_i: str | None = None) -> Recording:
    """TDMS files that already hold a conductance channel (e.g. the *_M.tdms
    files written by the EC4 software: G in siemens from i1 and E1-E2).
    The conductance is used as it is. A current channel, if there is one,
    is kept as the current (linear gain a); otherwise I = G G0 V."""
    g = np.asarray(td.channels[name_g], float)
    unit = cs.conductance_unit
    if unit == "auto":
        unit = str(td.channel_props[name_g].get("unit_string", "S")).strip()
        if unit not in CONDUCTANCE_UNITS + ["G0"]:
            unit = "S"
    G = np.clip(conductance_from_unit(g, unit), None, cs.g_ceiling)
    G = apply_filter(G, cs.conductance_filter)
    E, bias_name = _tdms_bias(td, cs, len(G))
    if name_i is not None:
        gain = cs.gain_a if cs.gain_a is not None else _scale(td.channel_props[name_i], 1.0)
        current = (td.channels[name_i] - cs.offset_a) * gain * cs.corr_a * cs.amplification
    else:
        current = G * G0 * E
    meta = {"source": td.path, "precomputed_G": True, "conductance_channel": name_g, "G_unit": unit,
            "bias_channel": bias_name,
            "file_properties": {k: str(v) for k, v in td.file_props.items()},
            "conversion": cs.to_dict()}
    return Recording(name=os.path.splitext(os.path.basename(td.path))[0], fs=td.fs, G=G,
                     piezo=_tdms_piezo(td, cs), bias=E, current=current, meta=meta, events=td.events)


def back_calculated_current(rec: Recording) -> np.ndarray:
    """Current that would flow without the series resistor: I_bc = G * G0 * E.

    This is column 1 of the .cvr format.
    """
    return rec.G * G0 * rec.bias_array()


# Column files (.cvr / .txt / .csv / Igor .pxp / .ibw) -> Recording

@dataclass
class ColumnSettings:
    """Which column of a text file, or which wave of an Igor file, holds what.

    A column is given by its 0-based index or by its name (header or wave
    name); empty = not in the file. Whatever the file gives is used as it
    is: a conductance column is not recalculated, a voltage column replaces
    the constant bias, a time column (or the Igor x scaling) sets the
    sampling rate. Only what is missing is calculated:
    G = I / (V - I R) / G0 from the current, I = G G0 V from the conductance.
    A '*' in the conductance or Ia column makes one recording per matching
    column (e.g. one Igor wave per trace); other columns with '*' are paired
    by the same text (G_*, Z_* -> G_12 with Z_12).
    """

    conductance_column: str = ""      # G given directly (empty = calculate it)
    conductance_unit: str = "G/G0"    # "G/G0" | "S" | "log G/G0"
    amplifier: str = "single stage"   # "single stage" (Ia only) | "dual stage" (Ia + Ib)
    current_a_column: str = "0"       # Ia, low gain, or the only current
    current_b_column: str = ""        # Ib, high gain (dual stage)
    current_scale_a: float = 1.0      # A per column unit (1 if the column is in A)
    current_scale_b: float = 1.0
    b_saturation: float | None = None  # |Ib| (column units) from which Ib is saturated; None = 99.5 % of max
    a_saturation: float | None = None  # |Ia| (column units) from which G = saturated G; None = off
    switch_logG: float | None = None   # optional: use Ib only below this log10(G/G0)
    voltage_column: str = ""          # measured bias (V); replaces the constant bias
    bias: float = 0.1                 # V, only used without a voltage column
    series_resistance: float = 0.0    # ohm, for G from the current
    g_ceiling: float = 100.0          # G/G0 for saturated samples
    piezo_column: str = "1"           # piezo / distance (empty = none; distance from sample index)
    piezo_scale: float = 1.0          # multiply the piezo column to get nm
    time_column: str = ""             # time (s) -> sampling rate
    fs: float = 10000.0               # Hz, if neither a time column nor the file gives it

    def to_dict(self) -> dict:
        return asdict(self)


def _wildcard(ref: str):
    return re.compile("^" + re.escape(ref.strip()).replace(r"\*", "(.*)") + "$", re.IGNORECASE)


def _column(table, ref: str, key: str | None, lenient: bool = False):
    """Array for a column reference; '*' is replaced by key. lenient: missing -> None."""
    ref = (ref or "").strip()
    if not ref:
        return None
    if "*" in ref:
        if key is None:
            raise ValueError(f"'{ref}' has a '*' but the conductance / Ia column has none.")
        ref = ref.replace("*", key)
        lenient = True
    try:
        return table.get(ref)
    except ValueError:
        if lenient:
            return None
        raise


def _one_recording(table, cs: ColumnSettings, name: str, key: str | None) -> Recording:
    g_col = _column(table, cs.conductance_column, key)
    ia_col = _column(table, cs.current_a_column, key)
    ib_col = _column(table, cs.current_b_column, key) if cs.amplifier == "dual stage" else None
    v_col = _column(table, cs.voltage_column, key)
    z_col = _column(table, cs.piezo_column, key, lenient=True)
    t_col = _column(table, cs.time_column, key)
    if g_col is None and ia_col is None and ib_col is None:
        raise ValueError("Give a conductance column or a current column.")
    # columns can differ in length (Igor waves): use the common part
    used = [c for c in (g_col, ia_col, ib_col, v_col, z_col, t_col) if c is not None]
    n = min(len(c) for c in used)
    g_col, ia_col, ib_col, v_col, z_col, t_col = (None if c is None else c[:n]
                                                  for c in (g_col, ia_col, ib_col, v_col, z_col, t_col))

    given, calculated = [], []
    bias = cs.bias
    if v_col is not None:
        bias = v_col
        given.append("voltage")
    ia = None if ia_col is None else ia_col * cs.current_scale_a
    ib = None if ib_col is None else ib_col * cs.current_scale_b
    if g_col is not None:
        G = np.clip(conductance_from_unit(g_col, cs.conductance_unit), None, cs.g_ceiling)
        given.append("conductance")
        if ia is not None or ib is not None:
            current = ia if ia is not None else ib
            given.append("current")
        else:
            current = G * G0 * bias
            calculated.append("current")
        frac_b = None
    else:
        b_ok = None
        if ib_col is not None:
            lim = cs.b_saturation
            if lim is None:
                lim = 0.995 * np.nanmax(np.abs(ib_col))
            b_ok = np.abs(ib_col) < lim
        a_sat = None if (cs.a_saturation is None or ia_col is None) else np.abs(ia_col) >= cs.a_saturation
        G, current, frac_b = currents_to_conductance(
            ia, ib, bias, cs.series_resistance, amplifier=cs.amplifier, b_ok=b_ok,
            switch_logG=cs.switch_logG, a_saturated=a_sat, g_ceiling=cs.g_ceiling)
        given.append("current")
        calculated.append("conductance")
    piezo = None
    if z_col is not None:
        piezo = z_col * cs.piezo_scale
        if np.ptp(piezo) == 0:
            piezo = None
        else:
            given.append("piezo")
    fs = table.fs or cs.fs
    if t_col is not None:
        dt = np.median(np.diff(t_col))
        if dt > 0:
            fs = 1.0 / dt
            given.append("time")
    G = np.abs(np.nan_to_num(G, nan=1e-12, posinf=cs.g_ceiling))
    return Recording(name=name, fs=fs, G=G, piezo=piezo, bias=bias, current=current,
                     meta={"columns": cs.to_dict(), "n_columns": len(table.columns),
                           "from_file": given, "calculated": calculated, "amplifier": cs.amplifier,
                           "fraction_high_gain": frac_b})


def table_to_recordings(table, cs: ColumnSettings, name: str = "data") -> list[Recording]:
    """Recordings from a :class:`bj_io.ColumnTable` (several if the main column has a '*')."""
    main = cs.conductance_column.strip() or cs.current_a_column.strip()
    if "*" not in main:
        return [_one_recording(table, cs, name, None)]
    pat = _wildcard(main)
    keys = [m.group(1) for m in (pat.match(c) for c in table.names) if m]
    if not keys:
        raise ValueError(f"No column matches '{main}'; available: {', '.join(table.names[:30])}")
    return [_one_recording(table, cs, f"{name}_{k}", k) for k in keys]


def columns_to_recording(arr: np.ndarray, cs: ColumnSettings, name: str = "data") -> Recording:
    """A plain 2-D array (columns 0, 1, ...) -> Recording."""
    from .bj_io import ColumnTable
    table = ColumnTable(name, {str(i): arr[:, i] for i in range(arr.shape[1])})
    return table_to_recordings(table, cs, name)[0]


def load_recordings(path: str, column_settings: ColumnSettings | None = None,
                    conversion_settings: ConversionSettings | None = None) -> list[Recording]:
    """Load any supported file; text and Igor files can give several recordings."""
    from . import bj_io

    ext = os.path.splitext(path)[1].lower()
    name = os.path.splitext(os.path.basename(path))[0]
    if ext == ".tdms":
        recs = [tdms_to_recording(bj_io.read_tdms(path), conversion_settings or ConversionSettings())]
    else:
        recs = table_to_recordings(bj_io.read_table(path), column_settings or ColumnSettings(), name)
    for r in recs:
        r.meta.setdefault("source", path)
    return recs


def load_recording(path: str, column_settings: ColumnSettings | None = None,
                   conversion_settings: ConversionSettings | None = None) -> Recording:
    """Load a file into one :class:`Recording` (the first, if it gives several)."""
    return load_recordings(path, column_settings, conversion_settings)[0]
