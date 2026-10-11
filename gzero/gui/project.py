"""What's loaded right now: recordings, traces and the trace settings."""

from __future__ import annotations

from PySide6.QtCore import QObject, Signal

from .. import streaming as stm, traces as tr
from ..conversion import ColumnSettings, ConversionSettings, Recording


class Project(QObject):
    recordings_changed = Signal()
    traces_changed = Signal()       # new segmentation -> every analysis is stale
    selection_changed = Signal()    # good/bad toggles or new cluster labels
    status = Signal(str)

    def __init__(self):
        super().__init__()
        self.recordings: list[Recording] = []
        self.traces: list[tr.Trace] = []
        self.holds: list[tr.Trace] = []
        self.segment_settings = tr.SegmentSettings()
        self.column_settings = ColumnSettings()
        self.conversion_settings = ConversionSettings()
        self.stream_settings = stm.StreamSettings()
        self.n_clusters = 0

    # data
    def add_recordings(self, recs: list[Recording], replace: bool = False):
        if replace:
            self.recordings = []
        names = {r.name for r in self.recordings}
        for r in recs:
            base, k = r.name, 2
            while r.name in names:      # keep names unique (traces refer to them)
                r.name = f"{base}_{k}"
                k += 1
            names.add(r.name)
            self.recordings.append(r)
        self.recordings_changed.emit()

    def remove_recording(self, name: str):
        self.recordings = [r for r in self.recordings if r.name != name]
        self.traces = [t for t in self.traces if t.rec_name != name]
        self.holds = [h for h in self.holds if h.rec_name != name]
        self.recordings_changed.emit()
        self.traces_changed.emit()

    def recording(self, name: str) -> Recording | None:
        return next((r for r in self.recordings if r.name == name), None)

    def full_recording(self, name: str, t_from: float | None = None, t_to: float | None = None) -> Recording | None:
        """The recording with all its signals. For a streamed (big) TDMS file the part
        t_from .. t_to (s, None = start / end) is read and converted now."""
        rec = self.recording(name)
        if rec is None or not rec.meta.get("streamed"):
            return rec
        fs = rec.meta["fs_full"]
        a = int((t_from or 0.0) * fs)
        b = int(t_to * fs) if t_to is not None else None
        part = stm.read_range(rec.meta["source"], ConversionSettings.from_dict(rec.meta["conversion"]), a, b)
        part.name = rec.name
        return part

    # traces
    def segment(self, ss: tr.SegmentSettings):
        self.segment_settings = ss
        traces, holds = [], []
        for rec in self.recordings:
            if rec.meta.get("streamed"):
                # big TDMS file: read it again chunk by chunk
                t, h = stm.stream_traces(rec, ss, start_index=len(traces))
            else:
                t, h = tr.extract_traces(rec, ss, start_index=len(traces))
            for i, hh in enumerate(h):
                hh.index = len(holds) + i
            traces += t
            holds += h
        self.traces, self.holds = traces, holds
        self.n_clusters = 0
        self.traces_changed.emit()

    def selected(self, cluster: int | None = None) -> list[tr.Trace]:
        out = [t for t in self.traces if t.selected]
        if cluster is not None and cluster >= 0:
            out = [t for t in out if t.cluster == cluster]
        return out

    def summary_text(self) -> str:
        s = tr.summary(self.traces)
        txt = (f"{s['total']} traces, {s['selected']} accepted ({s['percentage']:.1f} %), "
               f"{len(self.holds)} holds")
        if s["rejected_by"]:
            txt += "\nRejected: " +", ".join(f"{k}: {v}" for k, v in s["rejected_by"].items())
        cut = sum(r.meta.get("traces_cut_by_margin", 0) for r in self.recordings)
        if cut:
            txt += (f"\n{cut} traces were longer than the margin of a file read in pieces and may be cut "
                    "short: make 'margin' longer (Large TDMS files).")
        return txt
