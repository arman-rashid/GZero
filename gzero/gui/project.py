"""What's loaded right now: recordings, traces and the trace settings."""

from __future__ import annotations

from PySide6.QtCore import QObject, Signal

from .. import traces as tr
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

    # traces
    def segment(self, ss: tr.SegmentSettings):
        self.segment_settings = ss
        traces, holds = [], []
        for rec in self.recordings:
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
        return txt
