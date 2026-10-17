"""Blinking (fixed-gap I(t)) and change-point tabs."""

from __future__ import annotations

import os

import numpy as np
from PySide6.QtWidgets import QComboBox, QHBoxLayout, QLabel, QSpinBox, QWidget

from .. import blinking as bk, changepoints as cp, conversion as cv, export as ex, figures as figs
from .base import AnalysisTab
from .tabs_noise import _recording_box
from .widgets import PlotPanel, SettingsForm, button_row, error_box, make_button, scroll_panel

BLINK_LABELS = {
    "source": "segments", "smooth_ms": "smooth (ms)", "baseline_s": "baseline window (s)",
    "baseline_percentile": "baseline percentile", "up_threshold": "up threshold (decades)",
    "down_threshold": "down threshold (decades)", "min_on_ms": "min on (ms)", "min_off_ms": "min off (ms)",
    "g_min": "blink log G from", "g_max": "blink log G to", "bins_per_decade": "bins / decade",
    "time_bins": "time bins", "align_before_ms": "aligned: before (ms)", "align_after_ms": "aligned: after (ms)",
}
BLINK_CHOICES = {"source": ["recording", "holds"]}
BLINK_TIPS = {
    "source": "recording: the whole recording chosen above (a fixed-gap I(t) measurement); for a big TDMS file "
              "only the part between from and to.\nholds: the holds found in Data & Traces.",
    "baseline_s": "The baseline (tunnelling without a molecule) is a low percentile of log G in blocks of this "
                  "length, so it follows slow drift of the gap.",
    "baseline_percentile": "Must be below the fraction of time without a molecule.",
    "up_threshold": "A blink starts when log G rises this far above the baseline.",
    "down_threshold": "... and ends when it falls below this (hysteresis).",
    "min_off_ms": "Shorter gaps between two blinks are closed.",
    "g_min": "Blinks whose level is outside this window (metal contact, noise) are dropped.",
}

CP_LABELS = {
    "g_high": "window: log G below", "g_low": "window: log G above", "min_points": "min segment (points)",
    "f_threshold": "Chow F threshold", "max_plateau_slope": "max plateau slope (dec/nm)",
    "penalty": "PELT penalty factor", "max_points": "PELT: max points per trace",
    "noise_floor": "PELT: noise at least (decades)", "min_step": "PELT: merge levels closer than (dec)",
    "min_level_nm": "PELT: min level length (nm)", "bins_per_decade": "bins / decade",
}
CP_TIPS = {
    "g_high": "The analysis starts at the first point below this (after the gold contact).",
    "g_low": "... and stops at the first point below this (the noise floor).",
    "f_threshold": "Chow statistic of two lines against one at the best break (Hamill et al. 2024).",
    "max_plateau_slope": "The part before the rupture must be flatter than this to count as a molecular plateau; "
                         "a tunnelling decay falls by beta/ln10 (4-9 decades/nm).",
    "penalty": "Cost per change point = factor x 2 ln(n) (BIC-like). Larger: fewer steps.",
    "max_points": "Traces are block-averaged to at most this many points before PELT (speed).",
    "noise_floor": "PELT uses the noise of the trace, but at least this, so smooth ramps aren't cut into many "
                   "levels.",
}


class BlinkingTab(AnalysisTab):
    title = "Blinking"
    reset_on_traces = {"res": None, "segs": list}
    reset_on_recordings = {"res": None, "segs": list}

    def __init__(self, project, main):
        super().__init__(project, main)
        self.rec = QComboBox()
        self.form = SettingsForm("Blink detection", bk.BlinkSettings(), BLINK_LABELS, BLINK_CHOICES, BLINK_TIPS)
        self.forms["blinking"] = self.form
        self.result_label = QLabel("")
        self.result_label.setWordWrap(True)
        self.example = QSpinBox()
        self.example.setRange(0, 0)
        self.example.valueChanged.connect(lambda _: self._draw())
        row = QWidget()
        rl = QHBoxLayout(row)
        rl.setContentsMargins(0, 0, 0, 0)
        rl.addWidget(QLabel("show segment"))
        rl.addWidget(self.example)
        info = QLabel("Fixed-gap I(t): with the gap held still, a molecule jumping in gives a step up of the "
                      "current and its leaving a step down (Haiss et al. 2003, 2006). The blink levels give "
                      "the molecular conductance, the on and off times the lifetime and the formation rate.")
        info.setWordWrap(True)
        self.layout_.addWidget(scroll_panel(
            _recording_box(self, self.rec), self.form,
            button_row(make_button("Detect blinks", self.run),
                       make_button("Add simulated", self.add_simulated,
                                   "A 20 s fake I(t) recording: level -3.5, tau_on 20 ms, tau_off 100 ms.")),
            row, self.result_label, info, width=360))
        self.plot = PlotPanel()
        self.layout_.addWidget(self.plot, 1)
        self.res = None
        self.segs = []

    def add_simulated(self):
        k = sum(r.name.startswith("blink sim") for r in self.project.recordings)
        y, fs, _ = bk.simulate(seed=k)
        rec = cv.Recording(name=f"blink sim {k + 1}", fs=fs, G=10.0 ** y, meta={"simulated": True})
        self.project.add_recordings([rec])
        self.rec.setCurrentText(rec.name)

    def run(self):
        s = self.form_value("blinking")
        if s is None:
            return
        if s.source == "holds":
            segs = [(h.logG, h.fs) for h in self.project.holds]
            names = [f"hold {h.index}" for h in self.project.holds]
            if not segs:
                error_box(self, self.title, "No holds. Detect the traces in Data & Traces first (holds need a "
                                            "piezo signal), or use segments = recording.")
                return
            job = lambda: (bk.analyse(segs, s, names), segs)  # noqa: E731
        else:
            name = self.rec.currentText()
            if not name:
                error_box(self, self.title, "Load a recording first.")
                return
            t0, t1 = self.t_from.value(), self.t_to.value()

            def job():
                rec = self.project.full_recording(name, t0, t1)
                seg = [(np.log10(np.clip(rec.G, 1e-10, None)), rec.fs)]
                return bk.analyse(seg, s, [name]), seg

        def done(out):
            self.res, self.segs = out
            first = self.res.blinks[0].segment if self.res.blinks else 0
            self.example.blockSignals(True)
            self.example.setRange(0, max(0, len(self.segs) - 1))
            self.example.setValue(first)
            self.example.blockSignals(False)
            self.result_label.setText("<br>".join(self.res.summary_lines()))
            self._draw()

        self.background(job, on_done=done, busy="Detecting blinks ...")

    def panels(self):
        if self.res is None:
            return []
        k = self.example.value()
        example = None
        if 0 <= k < len(self.segs):
            lg, fs = self.segs[k]
            bl, ys, base, _ = bk.detect(lg, fs, self.res.settings, k)
            nm = self.res.segment_names[k] if k < len(self.res.segment_names) else f"segment {k}"
            example = (np.arange(len(lg)) / fs, lg, ys, base, bl, f"{nm}: {len(bl)} blinks")
        return figs.blinking(self.res, example)

    def _draw(self):
        self.plot.show_panels(self.panels(), 3)

    def has_results(self):
        return self.res is not None

    def export(self, folder, prefix):
        if self.res is None:
            return []
        r = self.res
        b = r.blinks
        f1 = os.path.join(folder, f"{prefix}_blinks.txt")
        ex.save_columns(f1, ["segment", "start_s", "duration_ms", "logG", "baseline_logG", "height_decades"],
                        [x.segment for x in b], [x.start / self.segs[x.segment][1] for x in b],
                        [1e3 * x.duration_s for x in b], [x.logG for x in b], [x.baseline for x in b],
                        [x.height for x in b])
        files = [f1]
        if len(r.off_times):
            f2 = os.path.join(folder, f"{prefix}_blink_off_times.txt")
            ex.save_columns(f2, ["off_time_ms"], 1e3 * r.off_times)
            files.append(f2)
        if r.hist:
            f3 = os.path.join(folder, f"{prefix}_blink_histogram.txt")
            ex.save_columns(f3, ["logG", "blinks"], *r.hist)
            files.append(f3)
        files.append(ex.save_summary(folder, [x.replace("<b>", "").replace("</b>", "") for x in r.summary_lines()]
                                     + ["settings: " + str(r.settings.to_dict())], f"{prefix}_blinking_summary.txt"))
        return files


class ChangePointTab(AnalysisTab):
    title = "Change points"
    reset_on_traces = {"res": None, "traces": list, "lh": None}

    def __init__(self, project, main):
        super().__init__(project, main)
        self.form = SettingsForm("Change points in breaking traces", cp.ChangePointSettings(), CP_LABELS,
                                 None, CP_TIPS)
        self.forms["changepoints"] = self.form
        self.result_label = QLabel("")
        self.result_label.setWordWrap(True)
        self.example = QSpinBox()
        self.example.setRange(0, 0)
        self.example.valueChanged.connect(lambda _: self._draw())
        row = QWidget()
        rl = QHBoxLayout(row)
        rl.setContentsMargins(0, 0, 0, 0)
        rl.addWidget(QLabel("show trace"))
        rl.addWidget(self.example)
        info = QLabel("Rupture: the Chow test finds where two straight lines describe log G(z) better than one; "
                      "the plateau conductance is the median before that break (Hamill et al. 2024). "
                      "Steps: PELT splits every trace into flat levels, so traces with two or more steps show "
                      "up (Killick et al. 2012).")
        info.setWordWrap(True)
        self.layout_.addWidget(scroll_panel(self.form, make_button("Find change points", self.run), row,
                                            self.result_label, info, width=360))
        self.plot = PlotPanel()
        self.layout_.addWidget(self.plot, 1)
        self.res = None
        self.traces = []
        self.lh = None

    def run(self):
        s = self.form_value("changepoints")
        if s is None:
            return
        traces = self.project.selected()
        if not traces:
            error_box(self, self.title, "No accepted traces. Detect the traces in Data & Traces first.")
            return

        def done(res):
            self.res, self.traces = res, traces
            self.lh = cp.level_histogram(res)
            first = next((r.trace for r in res.ruptures if r.significant), 0)
            self.example.blockSignals(True)
            self.example.setRange(0, max(0, len(traces) - 1))
            self.example.setValue(first)
            self.example.blockSignals(False)
            self.result_label.setText("<br>".join(res.summary_lines()))
            self._draw()

        self.background(cp.analyse, traces, s, on_done=done, busy="Finding change points ...")

    def panels(self):
        if self.res is None:
            return []
        k = self.example.value()
        example = None
        if 0 <= k < len(self.traces):
            t = self.traces[k]
            s = self.res.settings
            r = cp.rupture(t, s, k)
            f = cp.steps(t, s, k)
            pm = cp.plateaus(f, s) if f is not None else []
            lab = f"trace {t.index}: " + ("junction" if r is not None and r.significant else "blank")
            example = (t.z, t.logG, r, f, pm, lab)
        return figs.changepoints(self.res, example, self.lh)

    def _draw(self):
        self.plot.show_panels(self.panels(), 3)

    def has_results(self):
        return self.res is not None

    def export(self, folder, prefix):
        if self.res is None:
            return []
        r, s = self.res, self.res.settings
        R = r.ruptures
        f1 = os.path.join(folder, f"{prefix}_ruptures.txt")
        ex.save_columns(f1, ["trace", "z_break_nm", "chow_F", "plateau_logG", "length_nm", "slope_before",
                             "slope_after", "junction"],
                        [r.index[x.trace] for x in R], [x.z_break for x in R], [x.F for x in R],
                        [x.logG_plateau for x in R], [x.length for x in R], [x.slope_before for x in R],
                        [x.slope_after for x in R], [int(x.significant) for x in R])
        rows = [(r.index[f.trace], j, lv, ln, sl, int(p)) for f in r.steps
                for j, (lv, ln, sl, p) in enumerate(zip(f.levels, f.lengths, f.slopes, cp.plateaus(f, s)))]
        f2 = os.path.join(folder, f"{prefix}_step_levels.txt")
        ex.save_columns(f2, ["trace", "level", "logG", "length_nm", "slope", "plateau"],
                        *(list(c) for c in zip(*rows)) if rows else [[]] * 6)
        c, h, d = self.lh
        f3 = os.path.join(folder, f"{prefix}_step_level_histogram.txt")
        ex.save_columns(f3, ["logG", "summed_length_nm"], c, h)
        files = [f1, f2, f3]
        files.append(ex.save_summary(folder, [x.replace("<b>", "").replace("</b>", "") for x in r.summary_lines()]
                                     + ["settings: " + str(s.to_dict())], f"{prefix}_changepoints_summary.txt"))
        return files
