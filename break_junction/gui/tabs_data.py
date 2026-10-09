"""Convert, Data & Traces and Simulation tabs."""

from __future__ import annotations

import os

import numpy as np
from PySide6.QtCore import Qt
from PySide6.QtGui import QKeySequence, QShortcut
from PySide6.QtWidgets import (
    QCheckBox, QFileDialog, QGroupBox, QHBoxLayout, QLabel, QListWidget, QListWidgetItem,
    QPushButton, QSpinBox, QSplitter, QVBoxLayout, QWidget, QAbstractItemView,
)

from .. import bj_io, conversion as cv, figures as figs, simulation as sim, traces as tr
from ..analysis import centers
from .base import AnalysisTab
from .widgets import PlotPanel, SettingsForm, button_row, decimate_for_plot, error_box, make_button, scroll_panel

FILTER_KINDS = ["none", "savgol", "median", "moving average"]
FILTER_CHOICES = {"kind": FILTER_KINDS}
FILTER_LABELS = {"kind": "filter", "side_points": "side points", "poly_order": "polynomial order",
                 "in_log": "apply on log"}

CONVERSION_LABELS = {
    "module": "module (M_E#, M_i#a/b)",
    "amplifier": "current amplifier",
    "current_a_channel": "Ia channel", "current_b_channel": "Ib channel",
    "conductance_channel": "conductance channel", "conductance_unit": "conductance unit",
    "piezo_channel": "piezo channel",
    "calibration": "calibration",
    "gain_a": "gain Ia (A/V)", "gain_b": "gain Ib (A/V)",
    "offset_a": "offset Ia (V)", "offset_b": "offset Ib (V)",
    "corr_a": "correction Ia", "corr_b": "correction Ib",
    "amplification": "amplification",
    "invert_current": "invert current",
    "combine": "combine Ia / Ib",
    "b_saturation_V": "Ib saturation (V)",
    "a_saturation_V": "Ia saturation (V)",
    "bias_channel": "bias channel",
    "switch_logG": "switch to Ib below log G",
    "bias_source": "bias",
    "fixed_bias": "fixed bias (V)",
    "series_compensation": "series R compensation",
    "series_resistance": "series R (ohm)",
    "g_ceiling": "saturated G (G0)",
    "piezo_scale": "piezo scale (nm/V)",
    "bias_filter": "bias filter", "current_filter": "current filter",
    "conductance_filter": "conductance filter",
    "bias_filter.labels": FILTER_LABELS, "current_filter.labels": FILTER_LABELS,
    "conductance_filter.labels": FILTER_LABELS,
}
CONVERSION_CHOICES = {
    "amplifier": cv.AMPLIFIERS, "conductance_unit": ["auto"] + cv.CONDUCTANCE_UNITS,
    "calibration": ["linear", "diode"], "combine": ["auto", "a", "b"],
    "bias_source": ["measured", "fixed"],
    "bias_filter.choices": FILTER_CHOICES, "current_filter.choices": FILTER_CHOICES,
    "conductance_filter.choices": FILTER_CHOICES,
}
CONVERSION_TIPS = {
    "amplifier": "dual stage: low-gain Ia and high-gain Ib are measured together and merged\n"
                 "single stage: one current (the Ia channel) only, nothing to merge",
    "current_a_channel": "auto: M_i<module>a. Type a channel name to use another one; a channel that is "
                         "already in A needs gain 1.",
    "current_b_channel": "auto: M_i<module>b. Only used with a dual-stage amplifier.",
    "conductance_channel": "auto: use the channel 'G' only if there are no current channels.\n"
                           "Type a channel name to take the conductance directly from it (no calculation);\n"
                           "'none' = always calculate it from the current.",
    "conductance_unit": "Unit of the conductance channel; auto reads it from the channel (G0 or S).",
    "piezo_channel": "auto: U_Piezo. Type another channel name, or 'none'.",
    "gain_a": "Leave empty to use the scale factor stored with the channel in the TDMS file.",
    "gain_b": "Leave empty to use the scale factor stored with the channel in the TDMS file.",
    "calibration": "linear: I = (V - offset) * gain * correction (module 1)\n"
                   "diode: piecewise polynomial from the Calibration_M# group (logarithmic amplifier, module 2)",
    "combine": "Dual stage only. auto: Ib wherever it is not saturated, Ia elsewhere; a / b: one stage only.",
    "series_resistance": "Leave empty to use the current-limiting resistor stored in the Setup group.",
    "switch_logG": "Optional: use Ib only below this log10(G/G0).",
    "a_saturation_V": "Above this channel voltage the contact conducts more than the low-gain amplifier "
                      "can measure;\nsuch samples are set to the saturated G. Leave empty to switch off.",
    "bias_channel": "auto: use M_E1-E2 if present (as the acquisition software does), otherwise "
                    "M_E<module>.\nYou can also type a channel name.",
    "bias_source": "measured: the bias channel (an array) is used and the fixed bias is ignored;\n"
                   "fixed: a constant bias. Without a bias channel the fixed bias is used.",
    "g_ceiling": "Conductance assigned where the series resistor carries the whole bias (closed contact).",
}

COLUMN_LABELS = {
    "conductance_column": "conductance column", "conductance_unit": "conductance unit",
    "amplifier": "current amplifier",
    "current_a_column": "current Ia column", "current_b_column": "current Ib column",
    "current_scale_a": "Ia scale (A per unit)", "current_scale_b": "Ib scale (A per unit)",
    "b_saturation": "Ib saturated from |Ib|", "a_saturation": "Ia saturated from |Ia|",
    "switch_logG": "use Ib below log G",
    "voltage_column": "voltage column", "bias": "constant bias (V)",
    "series_resistance": "series R (ohm)", "g_ceiling": "saturated G (G0)",
    "piezo_column": "piezo / distance column", "piezo_scale": "piezo scale (to nm)",
    "time_column": "time column (s)", "fs": "sampling rate (Hz)",
}
COLUMN_CHOICES = {"conductance_unit": cv.CONDUCTANCE_UNITS, "amplifier": list(reversed(cv.AMPLIFIERS))}
_COL = "Column number (0 = first) or name (header / Igor wave name). Empty = not in the file."
COLUMN_TIPS = {
    "conductance_column": _COL + "\nIf given, G is taken from it and not calculated. "
                                 "A '*' (e.g. G_*) loads every matching column as its own recording.",
    "amplifier": "single stage: one current column (Ia).\n"
                 "dual stage: low-gain Ia and high-gain Ib, merged: Ib wherever it is not saturated.",
    "current_a_column": _COL + "\nThe current (single stage) or the low-gain stage Ia (dual stage).",
    "current_b_column": _COL + "\nHigh-gain stage Ib; only used with a dual-stage amplifier.",
    "current_scale_a": "Multiplies the Ia column to give A (1 if it is in A, the gain in A/V if it is in V).",
    "current_scale_b": "Multiplies the Ib column to give A.",
    "b_saturation": "|Ib| in column units from which Ib counts as saturated (Ia is used there).\n"
                    "Empty = 99.5 % of the largest |Ib| in the file.",
    "a_saturation": "|Ia| in column units from which the contact is beyond the range: G = saturated G. "
                    "Empty = off.",
    "voltage_column": _COL + "\nMeasured bias in V; replaces the constant bias.",
    "bias": "Used only without a voltage column. For .cvr files G = I_bc / V_bias / G0.",
    "series_resistance": "Series (current-limiting) resistor; G = I / (V - I R) / G0.",
    "piezo_column": _COL + "\nWithout piezo the distance comes from the point number x piezo rate.",
    "time_column": _COL + "\nSets the sampling rate (otherwise the Igor x scaling or the value below).",
    "fs": "Used if neither a time column nor the Igor x scaling gives the sampling rate.",
}


def link_amplifier(form, dual_only):
    """Grey out the high-gain (Ib) settings while the form is set to a single-stage amplifier."""
    combo = form.widgets["amplifier"]

    def update():
        for f in dual_only:
            form.widgets[f].setEnabled(combo.currentText() == "dual stage")

    def set_value(instance, _set=form.set_value):   # set_value blocks the combo's signals
        _set(instance)
        update()

    combo.currentTextChanged.connect(update)
    form.set_value = set_value
    update()

SEGMENT_LABELS = {
    "method": "segmentation", "use": "use traces", "smooth_ms": "piezo smoothing (ms)",
    "hold_fraction": "hold if speed below (fraction)", "min_segment_ms": "min segment (ms)",
    "classify": "breaking direction", "join_holds": "join pull-hold-pull",
    "distance_from": "distance from", "displacement_ratio": "displacement ratio (nm/nm)",
    "piezo_rate": "piezo rate (nm/point)",
    "zero_set": "zero set (log G)", "high_check": "high limit to check (log G)",
    "low_check": "low limit to check (log G)", "low_cut": "low limit to cut (log G)",
    "additional_points": "additional length (points)", "high_start": "start at (log G)",
    "noise_check": "reject noisy start", "noise_points": "noise points",
    "noise_limit": "noise limit (std log G)", "max_floor_noise": "max floor noise (std)",
}
SEGMENT_CHOICES = {"method": ["piezo", "threshold"], "use": ["break", "make", "both"],
                   "distance_from": ["piezo", "index"],
                   "classify": ["conductance", "piezo+", "piezo-"]}
SEGMENT_TIPS = {
    "method": "piezo: split where the piezo changes direction or stops.\n"
              "threshold: for files without a piezo signal; split at crossings of the high limit.",
    "classify": "conductance: a segment whose conductance falls is a breaking trace.\n"
                "piezo+ / piezo-: breaking traces are those where the piezo increases / decreases.",
    "join_holds": "Hold protocols stop the piezo in the middle of a pull. Join the pulls before and after "
                  "the hold into one trace; the hold samples are kept separately for the noise analysis.",
    "displacement_ratio": "MCBJ attenuation factor: junction displacement per nm of piezo elongation.",
    "zero_set": "Distance is set to zero where the trace first drops below this value.",
    "piezo_rate": "Distance per data point (nm), used for 'index' distance and for files "
                  "without a piezo signal.",
    "distance_from": "piezo: smoothed piezo signal x displacement ratio\n"
                     "index: point number x piezo rate",
}


# Convert tab

class ConvertTab(AnalysisTab):
    title = "Convert TDMS"

    def __init__(self, project, main):
        super().__init__(project, main)
        self.files = QListWidget()
        self.files.setSelectionMode(QAbstractItemView.ExtendedSelection)
        self.files.currentRowChanged.connect(self._show_info)
        self.form = SettingsForm("Conversion", project.conversion_settings,
                                 CONVERSION_LABELS, CONVERSION_CHOICES, CONVERSION_TIPS)
        self.forms["conversion"] = self.form
        link_amplifier(self.form, ("current_b_channel", "gain_b", "offset_b", "corr_b", "combine",
                                   "b_saturation_V", "switch_logG"))
        self.info = QLabel("Add .tdms files to convert.")
        self.info.setWordWrap(True)
        self.info.setTextInteractionFlags(Qt.TextSelectableByMouse)
        self.also_csv = QCheckBox("Also write .csv files (t, G, log G, piezo, bias, current)")

        files_box = QGroupBox("TDMS files")
        fl = QVBoxLayout(files_box)
        fl.addWidget(self.files)
        fl.addWidget(button_row(make_button("Add files", self.add_files),
                                make_button("Add folder", self.add_folder),
                                make_button("Remove", self.remove_files)))
        actions = QGroupBox("Run")
        al = QVBoxLayout(actions)
        al.addWidget(make_button("Preview selected", self.preview))
        al.addWidget(make_button("Convert all -> .cvr ...", self.convert_all,
                                 "Write .cvr files (back-calculated current, piezo nm)"))
        al.addWidget(self.also_csv)
        al.addWidget(make_button("Load all into analysis", self.load_into_project,
                                 "Convert in memory and add to the Data & Traces tab"))

        self.layout_.addWidget(scroll_panel(files_box, self.info, self.form, actions, width=360))
        self.plot = PlotPanel()
        self.layout_.addWidget(self.plot, 1)
        self._tdms_cache = {}

    # files
    def add_files(self):
        paths, _ = QFileDialog.getOpenFileNames(self, "TDMS files", self.main.last_dir, "TDMS (*.tdms)")
        self._add(paths)

    def add_folder(self):
        d = QFileDialog.getExistingDirectory(self, "Folder with TDMS files", self.main.last_dir)
        if d:
            self._add(bj_io.list_data_files(d, (".tdms",)))

    def _add(self, paths):
        existing = {self.files.item(i).data(Qt.UserRole) for i in range(self.files.count())}
        for p in paths:
            if p in existing:
                continue
            it = QListWidgetItem(os.path.basename(p))
            it.setData(Qt.UserRole, p)
            it.setToolTip(p)
            self.files.addItem(it)
        if paths:
            self.main.last_dir = os.path.dirname(paths[0])
            if self.files.currentRow() < 0:
                self.files.setCurrentRow(0)

    def remove_files(self):
        for it in self.files.selectedItems():
            self.files.takeItem(self.files.row(it))

    def paths(self):
        return [self.files.item(i).data(Qt.UserRole) for i in range(self.files.count())]

    def _tdms(self, path):
        if path not in self._tdms_cache:
            self._tdms_cache = {path: bj_io.read_tdms(path)}   # keep only one in memory
        return self._tdms_cache[path]

    def _show_info(self, row):
        if row < 0:
            return
        path = self.files.item(row).data(Qt.UserRole)
        try:
            td = self._tdms(path)
        except Exception as exc:  # noqa: BLE001
            self.info.setText(f"Cannot read {path}:\n{exc}")
            return
        mods = td.modules()
        lines = [
            f"<b>{os.path.basename(path)}</b>",
            f"{td.n_samples} samples at {td.fs:g} Hz ({td.n_samples / td.fs:.1f} s)",
            f"modules: {mods}; channels: {', '.join(td.channels)}",
            f"series resistor: {td.series_resistance():g} ohm",
            f"piezo modulation: {'on' if td.piezo_modulation() else 'off'}",
            f"events: {len(td.events)}",
        ]
        for m, cal in td.calibration.items():
            n = sum(1 for v in cal.values() if isinstance(v, list) and len(v) > 1)
            if n:
                lines.append(f"Calibration_M{m}: {n} polynomials (diode calibration available)")
        for name, props in td.channel_props.items():
            sf = props.get("scale_factor")
            if sf is not None:
                lines.append(f"{name}: scale {float(sf):g} {props.get('scale_factor_unit', '')}"
                             f" - {props.get('Description', '')}")
        self.info.setText("<br>".join(lines))
        if mods and self.form.widgets["module"].value() not in mods:
            self.form.widgets["module"].setValue(mods[0])

    # actions
    def _convert(self, path, cs):
        return cv.tdms_to_recording(bj_io.read_tdms(path), cs)

    def preview(self):
        row = self.files.currentRow()
        if row < 0:
            error_box(self, self.title, "Select a TDMS file first.")
            return
        cs = self.form_value("conversion")
        if cs is None:
            return
        path = self.files.item(row).data(Qt.UserRole)
        td = self._tdms(path)
        self.background(lambda: (td, cv.tdms_to_recording(td, cs)), on_done=self._plot_preview,
                        busy="Converting ...")

    def _plot_preview(self, res):
        td, rec = res
        fig = self.plot.figure
        fig.clear()
        ax = fig.subplots(3, 1, sharex=True)
        t = rec.time
        raw = rec.meta.get("current_channels", []) + [rec.meta.get("conductance_channel")]
        for name, color in zip(raw, ("tab:blue", "tab:orange", "tab:green")):
            if name in td.channels:
                x, y = decimate_for_plot(t, td.channels[name])
                ax[0].plot(x, y, lw=0.5, color=color, label=name)
        ax[0].set_ylabel("raw current\nchannels (V)")
        ax[0].legend(loc="upper right", fontsize=8)
        x, y = decimate_for_plot(t, rec.logG())
        ax[1].plot(x, y, lw=0.5, color="k")
        ax[1].set_ylabel("log(G/G$_0$)")
        for ts, lab in rec.events:
            c = {"Retract": "tab:red", "Approch": "tab:green", "Hold": "0.6"}.get(lab)
            if c and ts <= t[-1]:
                ax[1].axvline(ts, color=c, lw=0.6, alpha=0.6)
        if rec.piezo is not None:
            x, y = decimate_for_plot(t, rec.piezo)
            ax[2].plot(x, y, lw=0.7, color="tab:purple")
        ax[2].set_ylabel("piezo (nm)")
        ax[2].set_xlabel("time (s)")
        fhg = rec.meta.get("fraction_high_gain")
        fig.suptitle(f"{rec.name}  |  R = {rec.meta.get('series_resistance_ohm', 0):g} ohm"
                     f"  |  bias: {rec.meta.get('bias_channel', 'fixed')}"
                     + (f"  |  high-gain channel used {100 * fhg:.0f} %" if fhg is not None else "")
                     + "  |  event lines: green = approach, red = retract, grey = hold", fontsize=9)
        self.plot.draw()

    def convert_all(self):
        paths = self.paths()
        if not paths:
            error_box(self, self.title, "Add TDMS files first.")
            return
        cs = self.form_value("conversion")
        if cs is None:
            return
        out = QFileDialog.getExistingDirectory(self, "Output folder for .cvr files",
                                               os.path.dirname(paths[0]))
        if not out:
            return
        csv = self.also_csv.isChecked()

        def job():
            written = []
            for p in paths:
                rec = self._convert(p, cs)
                base = os.path.join(out, rec.name)
                bj_io.write_cvr(base + ".cvr", cv.back_calculated_current(rec),
                                rec.piezo if rec.piezo is not None else np.zeros(rec.n))
                written.append(base + ".cvr")
                if csv:
                    import pandas as pd
                    pd.DataFrame({"t_s": rec.time, "G_G0": rec.G, "logG": rec.logG(),
                                  "piezo_nm": rec.piezo if rec.piezo is not None else np.nan,
                                  "bias_V": rec.bias_array(), "current_A": rec.current}
                                 ).to_csv(base + ".csv", index=False, float_format="%.6g")
                    written.append(base + ".csv")
            return written

        self.background(job, on_done=lambda w: self.status(f"Wrote {len(w)} files to {out}"),
                        busy=f"Converting {len(paths)} files ...")

    def load_into_project(self):
        paths = self.paths()
        if not paths:
            error_box(self, self.title, "Add TDMS files first.")
            return
        cs = self.form_value("conversion")
        if cs is None:
            return
        self.project.conversion_settings = cs

        def job():
            return [self._convert(p, cs) for p in paths]

        def done(recs):
            self.project.add_recordings(recs)
            self.main.show_tab("Data & Traces")
            self.status(f"Loaded {len(recs)} recordings. Detect the traces in the Data & Traces tab.")

        self.background(job, on_done=done, busy="Converting ...")


# Data & Traces tab (loading, cutting, trace browser)

class DataTab(AnalysisTab):
    title = "Data & Traces"

    def __init__(self, project, main):
        super().__init__(project, main)
        self.rec_list = QListWidget()
        self.rec_list.setSelectionMode(QAbstractItemView.ExtendedSelection)
        self.rec_list.currentRowChanged.connect(lambda _: self.plot_overview())
        self.col_form = SettingsForm("Columns (text, .cvr, Igor .pxp / .ibw)", project.column_settings,
                                     COLUMN_LABELS, COLUMN_CHOICES, COLUMN_TIPS)
        self.seg_form = SettingsForm("Trace detection and cutting", project.segment_settings,
                                     SEGMENT_LABELS, SEGMENT_CHOICES, SEGMENT_TIPS)
        self.forms.update({"columns": self.col_form, "segmentation": self.seg_form})
        link_amplifier(self.col_form, ("current_b_column", "current_scale_b", "b_saturation", "switch_logG"))
        self.summary = QLabel("No data loaded.")
        self.summary.setWordWrap(True)

        data_box = QGroupBox("Recordings")
        dl = QVBoxLayout(data_box)
        dl.addWidget(self.rec_list)
        dl.addWidget(button_row(make_button("Load files", self.load_files),
                                make_button("Load folder", self.load_folder)))
        dl.addWidget(button_row(make_button("Remove", self.remove_selected),
                                make_button("Clear", self.clear_all),
                                make_button("Show columns", self.show_columns,
                                            "List the columns / waves of a text or Igor file")))
        self.seg_button = make_button("Detect / cut traces", self.segment)
        self.seg_button.setStyleSheet("font-weight: bold;")

        self.layout_.addWidget(scroll_panel(data_box, self.col_form, self.seg_form,
                                            self.seg_button, self.summary, width=350))

        # right side: overview + trace browser
        right = QSplitter(Qt.Vertical)
        self.overview = PlotPanel(figsize=(8, 3))
        right.addWidget(self.overview)

        browser = QWidget()
        bl = QVBoxLayout(browser)
        bl.setContentsMargins(0, 0, 0, 0)
        ctl = QHBoxLayout()
        self.idx = QSpinBox()
        self.idx.setRange(0, 0)
        self.idx.valueChanged.connect(self.plot_trace)
        self.only_sel = QCheckBox("Accepted only")
        self.only_sel.toggled.connect(self._browser_range)
        self.n_overlay = QSpinBox()
        self.n_overlay.setRange(1, 5000)
        self.n_overlay.setValue(50)
        self.trace_info = QLabel("")
        for w in (QLabel("trace"), self.idx,
                  make_button("<", lambda: self.idx.setValue(self.idx.value() - 1)),
                  make_button(">", lambda: self.idx.setValue(self.idx.value() + 1)),
                  make_button("Good (G)", lambda: self.mark(True)),
                  make_button("Bad (B)", lambda: self.mark(False)),
                  make_button("Auto", lambda: self.mark(None)),
                  self.only_sel, QLabel("overlay"), self.n_overlay,
                  make_button("Overlay", self.overlay)):
            ctl.addWidget(w)
        ctl.addStretch(1)
        bl.addLayout(ctl)
        bl.addWidget(self.trace_info)
        self.trace_plot = PlotPanel(figsize=(8, 3))
        bl.addWidget(self.trace_plot, 1)
        right.addWidget(browser)
        right.setSizes([350, 450])
        self.layout_.addWidget(right, 1)

        for key, fn in (("G", lambda: self.mark(True)), ("B", lambda: self.mark(False)),
                        ("Right", lambda: self.idx.setValue(self.idx.value() + 1)),
                        ("Left", lambda: self.idx.setValue(self.idx.value() - 1))):
            sc = QShortcut(QKeySequence(key), self)
            sc.setContext(Qt.WidgetWithChildrenShortcut)
            sc.activated.connect(fn)

        project.recordings_changed.connect(self._refresh_list)
        project.traces_changed.connect(self._after_segment)

    # loading
    def load_files(self):
        paths, _ = QFileDialog.getOpenFileNames(
            self, "Break-junction data", self.main.last_dir,
            "Data (*.cvr *.tdms *.txt *.csv *.dat *.tsv *.pxp *.ibw);;All files (*)")
        self._load(paths)

    def load_folder(self):
        d = QFileDialog.getExistingDirectory(self, "Folder with data files", self.main.last_dir)
        if d:
            self._load(bj_io.list_data_files(d))

    def _load(self, paths):
        if not paths:
            return
        self.main.last_dir = os.path.dirname(paths[0])
        cols = self.form_value("columns")
        if cols is None:
            return
        conv = self.main.tab("Convert TDMS").form_value("conversion")
        if conv is None:
            return
        self.project.column_settings = cols

        def job():
            recs, errors = [], []
            for p in paths:
                try:
                    recs.extend(cv.load_recordings(p, cols, conv))
                except Exception as exc:  # noqa: BLE001
                    errors.append(f"{os.path.basename(p)}: {exc}")
            return recs, errors

        def done(res):
            recs, errors = res
            self.project.add_recordings(recs)
            if errors:
                error_box(self, "Some files could not be read", "\n".join(errors))
            if recs and not self.project.traces:
                self.segment()

        self.background(job, on_done=done, busy=f"Loading {len(paths)} file(s) ...")

    def show_columns(self):
        path, _ = QFileDialog.getOpenFileName(
            self, "Show columns", self.main.last_dir,
            "Text / Igor (*.cvr *.txt *.csv *.dat *.tsv *.pxp *.ibw);;All files (*)")
        if not path:
            return

        def done(table):
            fs = f"; sampling rate in the file {table.fs:g} Hz" if table.fs else ""
            self.summary.setText(f"<b>{os.path.basename(path)}</b> ({len(table.columns)} columns{fs}):<br>"
                                 + table.describe(200))

        self.background(bj_io.read_table, path, on_done=done, busy="Reading ...")

    def remove_selected(self):
        for it in self.rec_list.selectedItems():
            self.project.remove_recording(it.data(Qt.UserRole))

    def clear_all(self):
        self.project.recordings, self.project.traces, self.project.holds = [], [], []
        self.project.recordings_changed.emit()
        self.project.traces_changed.emit()

    def _refresh_list(self):
        self.rec_list.blockSignals(True)
        self.rec_list.clear()
        for r in self.project.recordings:
            it = QListWidgetItem(f"{r.name}  ({r.n / r.fs:.1f} s, {'piezo' if r.piezo is not None else 'no piezo'})")
            it.setData(Qt.UserRole, r.name)
            tip = str(r.meta.get("source", ""))
            if r.meta.get("from_file"):
                tip += (f"\nfrom the file: {', '.join(r.meta['from_file'])}"
                        f"\ncalculated: {', '.join(r.meta['calculated']) or 'nothing'}")
            it.setToolTip(tip)
            self.rec_list.addItem(it)
        self.rec_list.blockSignals(False)
        if self.project.recordings:
            self.rec_list.setCurrentRow(0)
        self.plot_overview()

    # segment
    def segment(self):
        if not self.project.recordings:
            error_box(self, self.title, "Load data first.")
            return
        ss = self.form_value("segmentation")
        if ss is None:
            return
        self.background(lambda: self.project.segment(ss), busy="Detecting traces ...")

    def _after_segment(self):
        self.summary.setText(self.project.summary_text())
        self._browser_range()
        self.plot_overview()

    def _browser_list(self):
        return self.project.selected() if self.only_sel.isChecked() else self.project.traces

    def _browser_range(self):
        n = len(self._browser_list())
        self.idx.setRange(0, max(0, n - 1))
        self.plot_trace()

    # plots
    def plot_overview(self):
        fig = self.overview.figure
        fig.clear()
        item = self.rec_list.currentItem()
        rec = self.project.recording(item.data(Qt.UserRole)) if item else None
        if rec is None:
            self.overview.draw()
            return
        ax = fig.add_subplot(111)
        t = rec.time
        x, y = decimate_for_plot(t, rec.logG())
        ax.plot(x, y, lw=0.4, color="k")
        colors = {"break": "tab:red", "make": "tab:blue"}
        for tr_ in self.project.traces:
            if tr_.rec_name == rec.name:
                ax.axvspan(tr_.start / rec.fs, tr_.stop / rec.fs, color=colors.get(tr_.kind, "y"),
                           alpha=0.15 if tr_.selected else 0.05, lw=0)
        for h in self.project.holds:
            if h.rec_name == rec.name:
                ax.axvspan(h.start / rec.fs, h.stop / rec.fs, color="0.5", alpha=0.15, lw=0)
        ax.set_xlabel("time (s)")
        ax.set_ylabel("log(G/G$_0$)")
        if rec.piezo is not None:
            ax2 = ax.twinx()
            x, y = decimate_for_plot(t, rec.piezo)
            ax2.plot(x, y, lw=0.6, color="tab:purple", alpha=0.6)
            ax2.set_ylabel("piezo (nm)", color="tab:purple")
        ax.set_title(f"{rec.name}: red = breaking, blue = making, grey = hold (faint = rejected)",
                     fontsize=9)
        self.overview.draw()

    def current_trace(self):
        lst = self._browser_list()
        i = self.idx.value()
        return lst[i] if 0 <= i < len(lst) else None

    def plot_trace(self):
        fig = self.trace_plot.figure
        fig.clear()
        t = self.current_trace()
        if t is None:
            self.trace_info.setText("")
            self.trace_plot.draw()
            return
        ax, axh = fig.subplots(1, 2, sharey=True, gridspec_kw={"width_ratios": [4, 1]})
        color = "tab:green" if t.selected else "tab:red"
        ax.plot(t.z, t.logG, lw=0.7, color=color)
        ax.axvline(0, color="0.7", lw=0.5)
        ax.set_xlabel("distance (nm)")
        ax.set_ylabel("log(G/G$_0$)")
        edges = np.linspace(min(-7, t.logG.min()), max(1, t.logG.max()), 161)
        h = np.histogram(t.logG, edges)[0]
        axh.barh(centers(edges), h, height=edges[1] - edges[0], color="0.4")
        axh.set_xlabel("counts")
        state = "GOOD" if t.selected else f"BAD ({t.reason or 'manual'})"
        if t.manual is not None:
            state += " [manual]"
        self.trace_info.setText(f"trace {t.index} ({t.kind}) from {t.rec_name}: {len(t.logG)} points, "
                                f"{t.z[-1] - t.z[0]:.3f} nm  ->  {state}")
        self.trace_plot.draw()

    def mark(self, good):
        t = self.current_trace()
        if t is None:
            return
        t.manual = good
        self.summary.setText(self.project.summary_text())
        self.project.selection_changed.emit()
        if self.only_sel.isChecked():
            self._browser_range()
        else:
            self.plot_trace()
            self.idx.setValue(self.idx.value() + 1)

    def overlay(self):
        lst = self.project.selected()
        if not lst:
            return
        self.trace_info.setText("")
        self.trace_plot.show_panels(figs.trace_overlay(lst, self.n_overlay.value()), 1)

    def export(self, folder, prefix):
        from .. import export as ex
        files = []
        sel = self.project.selected()
        if sel:
            files.append(ex.save_traces(folder, sel, f"{prefix}_Good_Curves"))
            files.append(ex.save_traces_npz(folder, self.project.traces, f"{prefix}_all_traces"))
        return files

    def has_results(self):
        return bool(self.project.traces)

    def export_panels(self):
        return figs.trace_overlay(self.project.selected(), self.n_overlay.value())


# Simulation tab

SIM_LABELS = {
    "n_cycles": "number of traces", "fs": "sampling rate (Hz)", "speed": "speed (nm/s)",
    "displacement_ratio": "displacement ratio", "beta": "tunnelling decay beta (1/nm)",
    "contact_atoms": "contact size range (G0)", "atom_step": "atomic plateau length (nm)",
    "snap_back": "snap-back (nm)", "p_molecule": "molecule probability",
    "molecule_logG": "molecule log G", "molecule_sigma": "molecule spread (decades)",
    "plateau_length": "plateau length (nm)", "plateau_length_sd": "plateau length SD (nm)",
    "plateau_slope": "plateau slope (dec/nm)", "second_molecule_logG": "2nd molecule log G",
    "p_second": "2nd molecule fraction", "hold_s": "hold time (s)",
    "pull_distance": "pull before hold (nm)", "retract_distance": "retract after hold (nm)",
    "noise_exponent": "noise exponent n (NP ~ G^n)", "noise_amplitude": "flicker amplitude",
    "noise_floor_logG": "noise floor log G", "floor_noise": "floor noise",
    "measurement_noise": "white noise (relative)", "bias": "bias (V)", "seed": "random seed",
}


class SimulationTab(AnalysisTab):
    title = "Simulation"

    def __init__(self, project, main):
        super().__init__(project, main)
        self.form = SettingsForm("Synthetic MCBJ data", sim.SimulationSettings(), SIM_LABELS)
        self.forms["simulation"] = self.form
        self.replace = QCheckBox("Replace loaded data")
        self.info = QLabel("Generates a continuous recording (approach, pull, hold, retract) with a known "
                           "molecular conductance, plateau length and noise exponent, and then detects "
                           "the traces. Use it to check settings and analyses.")
        self.info.setWordWrap(True)
        self.layout_.addWidget(scroll_panel(self.info, self.form, self.replace,
                                            make_button("Generate", self.generate)))
        self.plot = PlotPanel()
        self.layout_.addWidget(self.plot, 1)

    def generate(self):
        s = self.form_value("simulation")
        if s is None:
            return
        replace = self.replace.isChecked()

        def done(rec):
            self.project.add_recordings([rec], replace=replace)
            self._plot(rec)
            data = self.main.tab("Data & Traces")
            data.segment()
            truth = rec.meta["truth"]
            mol = truth["molecule_logG"]
            self.status(f"Simulated {s.n_cycles} cycles. True mean molecular log G = "
                        f"{np.mean(mol):.2f}; true noise exponent n = {s.noise_exponent}.")

        self.background(sim.simulate, s, on_done=done, busy="Simulating ...")

    def _plot(self, rec):
        fig = self.plot.figure
        fig.clear()
        ax = fig.subplots(2, 1, sharex=True)
        n = min(rec.n, int(rec.fs * 30))
        t = rec.time[:n]
        x, y = decimate_for_plot(t, rec.logG()[:n])
        ax[0].plot(x, y, lw=0.5, color="k")
        ax[0].set_ylabel("log(G/G$_0$)")
        x, y = decimate_for_plot(t, rec.piezo[:n])
        ax[1].plot(x, y, lw=0.7, color="tab:purple")
        ax[1].set_ylabel("piezo (nm)")
        ax[1].set_xlabel("time (s)")
        ax[0].set_title("First 30 s of the simulated recording", fontsize=9)
        self.plot.draw()
