"""Reading and writing data files.

.tdms           raw files from the EC4 BJ software (BJRawData, Setup,
                Calibration_M#, Events groups)
.cvr            two columns, back-calculated current (A) and piezo (nm);
                some old files have a third column that isn't used
.txt/.csv/.dat  any numeric columns, the user says which one is which
.pxp / .ibw     Igor Pro packed experiment / binary wave (needs igor2); every
                numeric wave is a column, 2-D waves give one column per column
.lvm            LabVIEW measurement file (text); Delta_X in the header sets the
                sampling rate, the X_Value column becomes the column "time"
.xlsx / .xlsm   Excel workbook (openpyxl); every sheet is read, a header row
                gives the column names
.mat            MATLAB file; v4-v7.2 with scipy, v7.3 (HDF5) with h5py
.npy / .npz     NumPy array(s)
.h5 / .hdf5     HDF5 (h5py); every numeric data set is a column, a sampling
                rate is taken from an attribute fs / sampling_rate / dt
*.gz            any of the text formats compressed with gzip
Every format is turned into named numeric columns (ColumnTable), so the same
column settings work for all of them. 2-D arrays give one column per column
(the longer axis is taken as time).
"""

from __future__ import annotations

import os
import re
from dataclasses import dataclass, field

import numpy as np

TEXT_EXTENSIONS = (".cvr", ".txt", ".csv", ".dat", ".tsv", ".lvm")
GZIP_EXTENSIONS = tuple(e + ".gz" for e in TEXT_EXTENSIONS)
IGOR_EXTENSIONS = (".pxp", ".ibw")
EXCEL_EXTENSIONS = (".xlsx", ".xlsm")
MATLAB_EXTENSIONS = (".mat",)
NUMPY_EXTENSIONS = (".npy", ".npz")
HDF5_EXTENSIONS = (".h5", ".hdf5")
TABLE_EXTENSIONS = (TEXT_EXTENSIONS + GZIP_EXTENSIONS + IGOR_EXTENSIONS + EXCEL_EXTENSIONS + MATLAB_EXTENSIONS
                    + NUMPY_EXTENSIONS + HDF5_EXTENSIONS)
DATA_EXTENSIONS = (".tdms",) + TABLE_EXTENSIONS


def file_filter(tdms: bool = True) -> str:
    """Qt file-dialog filter with every readable format."""
    ext = DATA_EXTENSIONS if tdms else TABLE_EXTENSIONS
    return "Data (" + " ".join("*" + e for e in ext) + ");;All files (*)"


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
    n: int | None = None                             # samples (set when the channels are read lazily)

    @property
    def n_samples(self) -> int:
        if self.n is not None:
            return self.n
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


class LazyChannels(dict):
    """Channel name -> samples [start:stop], read from the file only when a channel is used."""

    def __init__(self, group, names, start, stop):
        super().__init__({n: None for n in names})
        self._group, self._start, self._stop = group, start, stop

    def __getitem__(self, name):
        v = dict.__getitem__(self, name)
        if v is None:
            v = np.asarray(self._group[name][self._start:self._stop], dtype=np.float64)
            dict.__setitem__(self, name, v)
        return v

    def get(self, name, default=None):
        return self[name] if name in self else default

    def values(self):
        return [self[k] for k in self]

    def items(self):
        return [(k, self[k]) for k in self]


class TdmsReader:
    """An EC4 BJ TDMS file opened for reading in pieces (needs npTDMS).

    Only the header and the small groups (Setup, Calibration, Events) are read
    when the file is opened; the raw channels are read slice by slice with
    :meth:`data`, so files much larger than the memory can be processed.

        with TdmsReader(path) as r:
            for start in range(0, r.n, 1_000_000):
                td = r.data(start, start + 1_000_000)     # a TdmsData of that range
    """

    def __init__(self, path: str):
        try:
            from nptdms import TdmsFile
        except ImportError as exc:  # pragma: no cover - environment dependent
            raise ImportError("Reading .tdms needs npTDMS:  pip install npTDMS") from exc
        self.path = path
        self._file = TdmsFile.open(path)
        f = self._file
        groups = {g.name: g for g in f.groups()}
        raw = groups.get("BJRawData")
        if raw is None:
            # fall back to the first group that holds numeric waveforms
            for g in f.groups():
                if any(np.dtype(c.dtype).kind == "f" for c in g.channels()):
                    raw = g
                    break
        if raw is None:
            self._file.close()
            raise ValueError(f"No raw data group found in {path}")
        self._raw = raw
        self.names, self.props, fs, lengths = [], {}, None, []
        for c in raw.channels():
            try:
                kind = np.dtype(c.dtype).kind
            except TypeError:
                continue
            if kind not in "fiu":
                continue
            self.names.append(c.name)
            self.props[c.name] = dict(c.properties)
            lengths.append(len(c))
            inc = c.properties.get("wf_increment")
            if fs is None and inc:
                fs = 1.0 / float(inc)
        self.n = int(min(lengths)) if lengths else 0
        self.fs = float(fs) if fs else 10000.0
        self.setup, self.calibration, self.events = _small_groups(groups)
        self.file_props = dict(f.properties)

    def data(self, start: int = 0, stop: int | None = None) -> TdmsData:
        stop = self.n if stop is None else min(int(stop), self.n)
        start = max(0, int(start))
        return TdmsData(path=self.path, channels=LazyChannels(self._raw, self.names, start, stop),
                        channel_props=self.props, fs=self.fs, setup=self.setup, calibration=self.calibration,
                        events=self.events, file_props=self.file_props, n=stop - start)

    def close(self):
        self._file.close()

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        self.close()


def _small_groups(groups):
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
    return setup, calibration, events


def tdms_info(path: str) -> TdmsData:
    """Header of a TDMS file without its samples (channels are read lazily if touched)."""
    r = TdmsReader(path)
    return r.data(0, r.n)


def read_tdms(path: str) -> TdmsData:
    """Read an EC4 BJ TDMS file completely (needs the npTDMS package)."""
    with TdmsReader(path) as r:
        td = r.data(0, r.n)
        td.channels = {k: v for k, v in td.channels.items()}   # read every channel now
        td.n = None
    return td


def _open_text(path: str):
    if path.lower().endswith(".gz"):
        import gzip
        return gzip.open(path, "rt", errors="replace")
    return open(path, "r", errors="replace")


def _sniff_delimiter(path: str, skip: int = 0):
    """Delimiter of the data lines (the first skip lines are the header)."""
    with _open_text(path) as fh:
        for k, line in enumerate(fh):
            if k < skip:
                continue
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
    with _open_text(path) as fh:
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

    skip = _count_header_lines(path)
    sep = _sniff_delimiter(path, skip)
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
    with _open_text(path) as fh:
        lines = [fh.readline() for _ in range(skip)]
    for line in reversed(lines):
        toks = [t.strip().strip('"') for t in re.split(r"[	,;]+|\s{2,}", line.strip()) if t.strip()]
        if len(toks) != ncol:
            toks = line.split()
        if len(toks) == ncol and len(set(toks)) == ncol:
            return toks
    return None


def _lvm_header(path: str):
    """(lines to skip, Delta_X, column names) of a LabVIEW .lvm file."""
    with _open_text(path) as fh:
        lines = [fh.readline() for _ in range(500)]
    ends = [i for i, line in enumerate(lines) if line.strip().startswith("***End_of_Header***")]
    if not ends:
        return 0, None, None
    skip, dx, names = ends[-1] + 1, None, None
    for line in lines[ends[0]:skip]:
        tok = line.rstrip("\r\n").split("\t")
        if tok[0] == "Delta_X":
            vals = [float(v) for v in tok[1:] if v.strip()]
            dx = vals[0] if vals else None
    hdr = lines[skip].rstrip("\r\n").split("\t") if skip < len(lines) else []
    if hdr and hdr[0].strip() in ("X_Value", "Time"):
        names = [h.strip() for h in hdr if h.strip() and h.strip() != "Comment"]
        skip += 1
    return skip, dx, names


def read_lvm(path: str) -> ColumnTable:
    """LabVIEW measurement file: tab separated, the header ends with ***End_of_Header***."""
    import pandas as pd

    skip, dx, names = _lvm_header(path)
    df = pd.read_csv(path, sep="\t", header=None, skiprows=skip, dtype=np.float64, on_bad_lines="skip",
                     engine="c")
    arr = df.to_numpy(np.float64)
    arr = arr[:, ~np.all(np.isnan(arr), axis=0)]
    if names is None or len(names) < arr.shape[1] or len(set(names)) != len(names):
        names = [str(i) for i in range(arr.shape[1])]
    cols = {n: arr[:, i] for i, n in enumerate(names[:arr.shape[1]])}
    fs = 1.0 / dx if dx and dx > 0 else None
    x = cols.pop("X_Value", None)
    if x is not None:
        steps = np.diff(x)
        if fs is None and len(x) > 1 and np.median(steps) > 0:
            fs = 1.0 / float(np.median(steps))
        cols["time"] = x
    return ColumnTable(path, cols, fs)


def _array_columns(name: str, a, out: dict):
    """Numeric array -> columns: 1-D as it is, 2-D one column per column (longer axis = time)."""
    a = np.asarray(a)
    if a.dtype.kind not in "fiubc" or a.size < 2:
        return
    if a.dtype.kind == "c":
        a = a.real
    a = np.squeeze(a.astype(np.float64))
    if a.ndim == 1:
        out[name] = a
    elif a.ndim == 2:
        if a.shape[0] < a.shape[1]:
            a = a.T
        for j in range(a.shape[1]):
            out[f"{name}[{j}]"] = a[:, j]


FS_KEYS = ("fs", "Fs", "sampling_rate", "SamplingRate", "sample_rate")
DT_KEYS = ("dt", "Delta_X", "wf_increment")


def _fs_from_attrs(attrs) -> float | None:
    """Sampling rate from attributes / scalar variables named fs, sampling_rate, dt, ..."""
    for keys, invert in ((FS_KEYS, False), (DT_KEYS, True)):
        for k in keys:
            if k in attrs:
                try:
                    v = float(np.ravel(attrs[k])[0])
                except (TypeError, ValueError, IndexError):
                    continue
                if v > 0:
                    return 1.0 / v if invert else v
    return None


def read_hdf5(path: str) -> ColumnTable:
    """Numeric data sets of an HDF5 file (also MATLAB v7.3); names are their paths (a:b:c)."""
    try:
        import h5py
    except ImportError as exc:  # pragma: no cover - environment dependent
        raise ImportError("Reading HDF5 / MATLAB v7.3 files needs h5py:  pip install h5py") from exc
    cols, scalars = {}, {}
    with h5py.File(path, "r") as f:
        fs = _fs_from_attrs(f.attrs)

        def visit(name, obj):
            if not isinstance(obj, h5py.Dataset) or obj.dtype.kind not in "fiu":
                return
            if obj.size == 1:
                scalars[name.split("/")[-1]] = obj[()]
            elif obj.ndim <= 2:
                _array_columns(name.replace("/", ":"), obj[()], cols)
                scalars.update({k: v for k, v in obj.attrs.items()})

        f.visititems(visit)
    if not cols:
        raise ValueError(f"No numeric data sets found in {os.path.basename(path)}")
    return ColumnTable(path, cols, fs or _fs_from_attrs(scalars))


def read_matlab(path: str) -> ColumnTable:
    """MATLAB .mat: every numeric variable (and struct field, as var.field) becomes column(s)."""
    from scipy.io import loadmat
    try:
        m = loadmat(path, squeeze_me=True, struct_as_record=False)
    except (NotImplementedError, ValueError):   # v7.3 files are HDF5
        with open(path, "rb") as fh:
            head = fh.read(1024)
        if b"\x89HDF\r\n\x1a\n" not in head:
            raise
        return read_hdf5(path)
    cols, scalars = {}, {}

    def walk(name, v):
        if hasattr(v, "_fieldnames"):
            for f in v._fieldnames:
                walk(f"{name}.{f}", getattr(v, f))
        elif np.ndim(v) == 0 and isinstance(v, (int, float, np.number)):
            scalars[name.split(".")[-1]] = v
        elif isinstance(v, np.ndarray):
            _array_columns(name, v, cols)

    for k, v in m.items():
        if not k.startswith("__"):
            walk(k, v)
    if not cols:
        raise ValueError(f"No numeric variables found in {os.path.basename(path)}")
    return ColumnTable(path, cols, _fs_from_attrs(scalars))


def read_numpy(path: str) -> ColumnTable:
    """.npy: one array (2-D gives the columns 0, 1, ...); .npz: every array by its name."""
    cols, scalars = {}, {}
    if path.lower().endswith(".npz"):
        with np.load(path, allow_pickle=False) as z:
            for k in z.files:
                a = z[k]
                if a.size == 1:
                    scalars[k] = a
                else:
                    _array_columns(k, a, cols)
    else:
        tmp = {}
        _array_columns("a", np.load(path, allow_pickle=False), tmp)
        cols = {("0" if k == "a" else k[2:-1]): v for k, v in tmp.items()}
    if not cols:
        raise ValueError(f"No numeric arrays found in {os.path.basename(path)}")
    return ColumnTable(path, cols, _fs_from_attrs(scalars))


def read_excel(path: str) -> ColumnTable:
    """Every sheet of an Excel workbook. The row above the first all-numeric row gives
    the column names; columns of the second and later sheets are named sheet:column."""
    import pandas as pd

    try:
        sheets = pd.read_excel(path, sheet_name=None, header=None)
    except ImportError as exc:  # pragma: no cover - environment dependent
        raise ImportError("Reading Excel files needs openpyxl:  pip install openpyxl") from exc
    cols = {}
    for k, (sheet, df) in enumerate(sheets.items()):
        df = df.dropna(how="all").dropna(axis=1, how="all")
        if df.empty:
            continue
        num = df.apply(pd.to_numeric, errors="coerce")
        full = num.notna().all(axis=1).to_numpy()
        first = int(np.argmax(full)) if full.any() else 0
        names = None
        if first > 0:
            hdr = [str(x).strip() for x in df.iloc[first - 1].tolist()]
            if len(set(hdr)) == len(hdr) and all(h and h != "nan" for h in hdr):
                names = hdr
        data = num.iloc[first:].to_numpy(np.float64)
        if names is None:
            names = [str(i) for i in range(data.shape[1])]
        pre = "" if k == 0 else f"{sheet}:"
        for j, n in enumerate(names):
            ok = np.flatnonzero(np.isfinite(data[:, j]))
            if len(ok) >= 2:
                cols[pre + n] = data[:ok[-1] + 1, j]
    if not cols:
        raise ValueError(f"No numeric columns found in {os.path.basename(path)}")
    return ColumnTable(path, cols)


def read_text_table(path: str) -> ColumnTable:
    if path.lower().removesuffix(".gz").endswith(".lvm"):
        return read_lvm(path)
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
    low = path.lower()
    if low.endswith(IGOR_EXTENSIONS):
        return read_igor(path)
    if low.endswith(EXCEL_EXTENSIONS):
        return read_excel(path)
    if low.endswith(MATLAB_EXTENSIONS):
        return read_matlab(path)
    if low.endswith(NUMPY_EXTENSIONS):
        return read_numpy(path)
    if low.endswith(HDF5_EXTENSIONS):
        return read_hdf5(path)
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
