"""Reading and writing data files.

.tdms           raw files from the EC4 BJ software (BJRawData, Setup,
                Calibration_M#, Events groups)
.cvr            two columns, back-calculated current (A) and piezo (nm);
                some old files have a third column that isn't used
.txt/.csv/.dat  any numeric columns, the user says which one is which
.pxp / .ibw     Igor Pro packed experiment / binary wave (needs igor2); every
                numeric wave is a column, 2-D waves give one column per column
"""

from __future__ import annotations

import os
import re
from dataclasses import dataclass, field

import numpy as np

TEXT_EXTENSIONS = (".cvr", ".txt", ".csv", ".dat", ".tsv")
IGOR_EXTENSIONS = (".pxp", ".ibw")
DATA_EXTENSIONS = (".tdms",) + TEXT_EXTENSIONS + IGOR_EXTENSIONS


@dataclass
class ColumnTable:
    """Named numeric columns from a text or Igor file (lengths may differ).

    fs is the sampling rate stored in the file (Igor x scaling), None if unknown.
    """

    path: str
    columns: dict           # name -> 1-D float64 array, in file order
    fs: float | None = None
    units: dict = field(default_factory=dict)

    @property
    def names(self) -> list[str]:
        return list(self.columns)

    def find(self, ref: str) -> str | None:
        """Column name for ref: an exact name, a case-insensitive name or a 0-based index."""
        ref = str(ref).strip()
        if not ref:
            return None
        if ref in self.columns:
            return ref
        low = {k.lower(): k for k in self.columns}
        if ref.lower() in low:
            return low[ref.lower()]
        if re.fullmatch(r"\d+", ref):
            i = int(ref)
            if i < len(self.columns):
                return self.names[i]
            raise ValueError(f"Column {i} is not in {os.path.basename(self.path)} "
                             f"(it has {len(self.columns)} columns, numbered from 0).")
        raise ValueError(f"No column '{ref}' in {os.path.basename(self.path)}; "
                         f"available: {', '.join(self.names[:30])}")

    def get(self, ref: str):
        name = self.find(ref)
        return None if name is None else self.columns[name]

    def describe(self, n: int = 40) -> str:
        parts = [f"{i}: {k} ({len(v)})" for i, (k, v) in enumerate(self.columns.items())]
        more = f" ... (+{len(parts) - n})" if len(parts) > n else ""
        return ", ".join(parts[:n]) + more


@dataclass
class TdmsData:
    """Everything useful from one EC4 BJ .tdms file."""

    path: str
    channels: dict          # name -> float64 array (raw volts as stored)
    channel_props: dict     # name -> dict of TDMS properties
    fs: float               # sampling rate (Hz)
    setup: dict = field(default_factory=dict)        # Setup group, item -> value string
    calibration: dict = field(default_factory=dict)  # module -> {item: [floats] or str}
    events: list = field(default_factory=list)       # [(t_seconds_from_first, label)]
    file_props: dict = field(default_factory=dict)

    @property
    def n_samples(self) -> int:
        return len(next(iter(self.channels.values()))) if self.channels else 0

    def modules(self) -> list[int]:
        """Measurement modules present (1, 2, ...) based on channel names."""
        mods = set()
        for name in self.channels:
            m = re.match(r"M_(?:E|i)(\d)", name)
            if m:
                mods.add(int(m.group(1)))
        return sorted(mods)

    def series_resistance(self) -> float:
        """Current-limiting (series) resistor in ohm, 0 if disabled/unknown."""
        state = self.setup.get("CurrentLimingResistor.State", "").strip().upper()
        value = self.setup.get("CurrentLimingResistor.Value", "")
        if state and state != "TRUE":
            return 0.0
        m = re.search(r"[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?", value)
        return float(m.group()) if m else 0.0

    def piezo_modulation(self) -> bool:
        return self.setup.get("Piezo_Modulation.Type", "Off").strip().lower() not in ("off", "")


def _parse_number_list(text: str):
    text = (text or "").strip()
    if not text:
        return []
    out = []
    for tok in text.replace("\r", "").replace("\n", "").split(";"):
        tok = tok.strip()
        if tok:
            try:
                out.append(float(tok))
            except ValueError:
                return text
    return out


def read_tdms(path: str) -> TdmsData:
    """Read an EC4 BJ TDMS file (needs the npTDMS package)."""
    try:
        from nptdms import TdmsFile
    except ImportError as exc:  # pragma: no cover - environment dependent
        raise ImportError("Reading .tdms needs npTDMS:  pip install npTDMS") from exc

    f = TdmsFile.read(path)
    groups = {g.name: g for g in f.groups()}
    raw = groups.get("BJRawData")
    if raw is None:
        # fall back to the first group that holds numeric waveforms
        for g in f.groups():
            if any(c.data.dtype.kind == "f" for c in g.channels()):
                raw = g
                break
    if raw is None:
        raise ValueError(f"No raw data group found in {path}")

    channels, props, fs = {}, {}, None
    for c in raw.channels():
        data = c[:]
        if data.dtype.kind not in "fiu":
            continue
        channels[c.name] = np.asarray(data, dtype=np.float64)
        props[c.name] = dict(c.properties)
        inc = c.properties.get("wf_increment")
        if fs is None and inc:
            fs = 1.0 / float(inc)

    setup = {}
    if "Setup" in groups:
        cols = [c[:] for c in groups["Setup"].channels()]
        # stored as a 2-row table: row 0 = item names, row 1 = values
        if cols and len(cols[0]) >= 2:
            for col in cols:
                setup[str(col[0])] = str(col[1])

    calibration = {}
    for gname, g in groups.items():
        m = re.match(r"Calibration_M(\d)", gname)
        if not m:
            continue
        names = {c.name: c[:] for c in g.channels()}
        if "Item" in names and "Value" in names:
            calibration[int(m.group(1))] = {
                str(k): _parse_number_list(str(v)) for k, v in zip(names["Item"], names["Value"])
            }

    events = []
    if "Events" in groups:
        ev = {c.name: c[:] for c in groups["Events"].channels()}
        if "EventTime" in ev:
            times = ev["EventTime"]
            labels = ev.get("EventType_Info", [""] * len(times))
            t0 = times[0]
            for t, lab in zip(times, labels):
                events.append(((t - t0) / np.timedelta64(1, "s"), str(lab)))

    return TdmsData(
        path=path,
        channels=channels,
        channel_props=props,
        fs=float(fs) if fs else 10000.0,
        setup=setup,
        calibration=calibration,
        events=events,
        file_props=dict(f.properties),
    )


def _sniff_delimiter(path: str):
    with open(path, "r", errors="replace") as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if "\t" in line:
                return "\t"
            if "," in line:
                return ","
            if ";" in line:
                return ";"
            return r"\s+"
    return r"\s+"


def _count_header_lines(path: str, max_check: int = 50) -> int:
    n = 0
    with open(path, "r", errors="replace") as fh:
        for line in fh:
            s = line.strip()
            if not s:
                n += 1
                continue
            first = re.split(r"[\t,;\s]+", s)[0]
            try:
                float(first)
                return n
            except ValueError:
                n += 1
            if n >= max_check:
                break
    return n


def read_columns(path: str) -> np.ndarray:
    """Fast numeric column reader (.cvr/.txt/.csv). Returns a 2-D float array."""
    import pandas as pd

    sep = _sniff_delimiter(path)
    skip = _count_header_lines(path)
    df = pd.read_csv(path, sep=sep, header=None, skiprows=skip, engine="c",
                     dtype=np.float64, on_bad_lines="skip")
    arr = df.to_numpy(dtype=np.float64)
    if arr.ndim == 1:
        arr = arr[:, None]
    # drop all-NaN trailing columns (lines ending in a delimiter)
    keep = ~np.all(np.isnan(arr), axis=0)
    return arr[:, keep]


def _header_names(path: str, skip: int, ncol: int):
    """Column names from the last header line, if it has one name per column."""
    if skip == 0:
        return None
    with open(path, "r", errors="replace") as fh:
        lines = [fh.readline() for _ in range(skip)]
    for line in reversed(lines):
        toks = [t.strip().strip('"') for t in re.split(r"[	,;]+|\s{2,}", line.strip()) if t.strip()]
        if len(toks) != ncol:
            toks = line.split()
        if len(toks) == ncol and len(set(toks)) == ncol:
            return toks
    return None


def read_text_table(path: str) -> ColumnTable:
    arr = read_columns(path)
    names = _header_names(path, _count_header_lines(path), arr.shape[1])
    if names is None:
        names = [str(i) for i in range(arr.shape[1])]
    return ColumnTable(path, {n: arr[:, i] for i, n in enumerate(names)})


def _igor_text(x) -> str:
    """Igor header text (bytes or an array of single characters) -> str."""
    if x is None:
        return ""
    if not isinstance(x, bytes):
        x = b"".join(np.asarray(x, dtype="S1").ravel().tolist())
    return x.split(b"\0")[0].decode("latin-1").strip()


def _igor_wave_columns(name: str, wave: dict, out: dict, units: dict):
    """Add the numeric columns of one Igor wave to out; return (x step, x unit)."""
    w = wave["wave"]
    data = w.get("wData")
    if data is None or not hasattr(data, "dtype") or data.dtype.kind not in "fiuc" or data.size < 2:
        return None, ""
    data = np.asarray(data)
    if data.dtype.kind == "c":
        data = data.real
    data = data.astype(np.float64)
    hdr = w.get("wave_header", {})
    try:
        if "sfA" in hdr:        # version 5
            dx, xunit = float(hdr["sfA"][0]), _igor_text(np.asarray(hdr["dimUnits"])[0])
        else:                   # versions 1-3
            dx, xunit = float(hdr.get("hsA", 1.0)), _igor_text(hdr.get("xUnits"))
        xunit = _igor_text(w.get("dimension_units")) or xunit
        unit = _igor_text(w.get("data_units")) or _igor_text(hdr.get("dataUnits"))
    except (TypeError, ValueError, IndexError):
        dx, xunit, unit = None, "", ""
    if data.ndim == 1:
        out[name] = data
        units[name] = unit
    else:
        data = data.reshape(data.shape[0], -1)
        for j in range(data.shape[1]):
            out[f"{name}[{j}]"] = data[:, j]
            units[f"{name}[{j}]"] = unit
    return dx, xunit


def read_igor(path: str) -> ColumnTable:
    """Read the numeric waves of an Igor .pxp or .ibw file (needs igor2).

    Waves inside data folders are named folder:wave; a 2-D wave gives the
    columns name[0], name[1], ... The sampling rate comes from the x scaling
    of the first wave with x unit "s" (or, if none has a unit, the first
    wave whose x step is not 1).
    """
    try:
        import igor2.binarywave
        import igor2.packed
    except ImportError as exc:  # pragma: no cover - environment dependent
        raise ImportError("Reading Igor files needs igor2:  pip install igor2") from exc

    cols, units, steps = {}, {}, []
    if path.lower().endswith(".ibw"):
        wave = igor2.binarywave.load(path)
        name = _igor_text(wave["wave"]["wave_header"].get("bname")) or "wave"
        steps.append(_igor_wave_columns(name, wave, cols, units))
    else:
        _, fs_tree = igor2.packed.load(path)

        def walk(folder: dict, prefix: str):
            for key, item in folder.items():
                k = key.decode("latin-1") if isinstance(key, bytes) else str(key)
                if isinstance(item, dict):
                    walk(item, f"{prefix}{k}:")
                elif hasattr(item, "wave"):
                    steps.append(_igor_wave_columns(prefix + k, item.wave, cols, units))

        walk(fs_tree.get("root", fs_tree), "")
    if not cols:
        raise ValueError(f"No numeric waves found in {os.path.basename(path)}")
    steps = [(dx, u) for dx, u in steps if dx and dx > 0]
    dx = next((dx for dx, u in steps if u == "s"), None)
    if dx is None and not any(u for _, u in steps):
        dx = next((dx for dx, _ in steps if dx != 1.0), None)
    return ColumnTable(path, cols, 1.0 / dx if dx else None, units)


def read_table(path: str) -> ColumnTable:
    """Any non-TDMS data file as a :class:`ColumnTable`."""
    if path.lower().endswith(IGOR_EXTENSIONS):
        return read_igor(path)
    return read_text_table(path)


def write_cvr(path: str, current_bc: np.ndarray, piezo_nm: np.ndarray) -> None:
    """Write a .cvr file (back-calculated current A, piezo nm)."""
    import pandas as pd

    df = pd.DataFrame({"i": np.asarray(current_bc, float), "z": np.asarray(piezo_nm, float)})
    df.to_csv(path, sep="\t", header=False, index=False, float_format="%.5E", lineterminator="\r\n")


def list_data_files(folder: str, extensions=DATA_EXTENSIONS) -> list[str]:
    files = [os.path.join(folder, f) for f in os.listdir(folder)
             if f.lower().endswith(tuple(extensions)) and os.path.isfile(os.path.join(folder, f))]
    return sorted(files)
