"""Thermopower (Seebeck) tab."""

from __future__ import annotations

import os

import numpy as np
from PySide6.QtWidgets import QComboBox, QFileDialog, QFormLayout, QGroupBox, QInputDialog, QLabel, QTableWidget, \
    QTableWidgetItem

from .. import bj_io, export as ex, figures as figs, style, thermo as th
from .base import AnalysisTab
from .widgets import PlotPanel, SettingsForm, button_row, error_box, make_button, scroll_panel

THERMO_LABELS = {
    "v_column": "voltage column", "voltage_unit": "voltage unit", "amplifier_gain": "amplifier gain (divide V by)",
    "invert_voltage": "invert voltage sign", "t_column": "temperature column", "t_mode": "temperature column is",
    "t_reference": "T reference (cold side)", "g_column": "conductance column (optional)",
    "g_unit": "conductance unit", "g_min": "keep log G from", "g_max": "keep log G to",
    "unit": "histogram of", "dT_step": "dT step (K)", "min_count": "min. samples / files per dT",
    "v_bins": "voltage bins", "s_electrode": "S of the electrodes (uV/K)", "s_bins": "Seebeck bins",
}
THERMO_CHOICES = {"voltage_unit": list(th.V_UNITS), "t_mode": th.T_MODES, "g_unit": th.G_UNITS,
                  "unit": ["samples", "file means"]}
THERMO_TIPS = {
    "v_column": "Column number (0 = first) or header name of the thermoelectric voltage.",
    "amplifier_gain": "If the voltage was recorded after an amplifier, its gain: V is divided by it.",
    "invert_voltage": "Use if your wiring measures V with the opposite sign (S comes out with the wrong sign "
                      "for a known molecule).",
    "t_mode": "temperature difference dT: the column holds dT.\n"
              "hot-side temperature: the column holds the heated electrode's temperature, dT = T - T reference.\n"
              "no column: you type dT for every data set when you load it.",
    "g_column": "Optional. With it, samples outside the log G window are left out and S is plotted against log G.",
    "unit": "samples: every point is one count (as in XMe). file means: one value per file and dT "
            "(one per junction, as in Reddy et al. 2007).",
    "dT_step": "dT values are grouped to multiples of this.",
    "s_electrode": "S = -dV/d(dT) + S_electrode. Leave 0 for S relative to the electrode leads, or enter the "
                   "Seebeck coefficient of your electrode metal.",
}


class ThermoTab(AnalysisTab):
    title = "Thermopower"

    def __init__(self, project, main):
        super().__init__(project, main)
        self.form = SettingsForm("Thermovoltage and temperature", th.ThermoSettings(), THERMO_LABELS, THERMO_CHOICES,
                                 THERMO_TIPS)
        self.forms["thermo"] = self.form
        self.datasets: list[th.ThermoDataset] = []
        self.results: dict = {}
        self.table = QTableWidget(0, 3)
        self.table.setHorizontalHeaderLabels(["data set", "files", "dT (K)"])
        self.table.setMinimumHeight(120)
        self.show_combo = QComboBox()
        self.show_combo.currentIndexChanged.connect(lambda _: self._draw())
        box = QGroupBox("Data sets")
        lay = QFormLayout(box)
        lay.addRow(button_row(make_button("Add files ...", self.add_files), make_button("Add folder ...", self.add_folder)))
        lay.addRow(button_row(make_button("Add simulated", self.add_simulated,
                                          "Fake traces at dT = 0, 10, 20, 30 K with S = 10 uV/K."),
                              make_button("Remove", self.remove)))
        lay.addRow(self.table)
        lay.addRow("show", self.show_combo)
        self.result_label = QLabel("Load files with the thermoelectric voltage and the temperature (difference). "
                                   "S = -dV/d(dT) from the most probable voltage at every dT.")
        self.result_label.setWordWrap(True)
        self.layout_.addWidget(scroll_panel(box, self.form, make_button("Analyse", self.run), self.result_label,
                                            width=380))
        self.plot = PlotPanel()
        self.layout_.addWidget(self.plot, 1)

    # data sets
    def _load(self, paths, name):
        s = self.form_value("thermo")
        if s is None or not paths:
            return
        dT = None
        if s.t_mode.startswith("no column"):
            dT, ok = QInputDialog.getDouble(self, "Temperature difference", f"dT for '{name}' (K):", 10.0, -1000,
                                            1000, 2)
            if not ok:
                return
        files, bad = [], []
        for p in paths:
            try:
                files.append(th.load_file(p, s, dT))
            except Exception as exc:  # noqa: BLE001
                bad.append(f"{os.path.basename(p)}: {exc}")
        if bad:
            error_box(self, self.title, "Some files could not be read:\n" + "\n".join(bad[:20]))
        if files:
            self.datasets.append(th.ThermoDataset(name, files))
            self._refresh()

    def add_files(self):
        paths, _ = QFileDialog.getOpenFileNames(self, "Thermovoltage files", self.main.last_dir,
                                                "Data (*.txt *.csv *.dat *.tsv *.pxp *.ibw);;All files (*)")
        if not paths:
            return
        self.main.last_dir = os.path.dirname(paths[0])
        name, ok = QInputDialog.getText(self, "Data set", "Name:", text=os.path.basename(os.path.dirname(paths[0])))
        if ok:
            self._load(paths, name.strip() or f"set {len(self.datasets) + 1}")

    def add_folder(self):
        folder = QFileDialog.getExistingDirectory(self, "Folder with thermovoltage files", self.main.last_dir)
        if not folder:
            return
        self.main.last_dir = folder
        paths = [p for p in bj_io.list_data_files(folder) if not p.lower().endswith(".tdms")]
        self._load(paths, os.path.basename(folder))

    def add_simulated(self):
        k = len(self.datasets)
        self.datasets.append(th.ThermoDataset(f"sim S=10 ({k + 1})", th.simulate(S=10.0, seed=k)))
        self._refresh()

    def remove(self):
        for r in sorted({i.row() for i in self.table.selectedIndexes()}, reverse=True):
            self.results.pop(self.datasets[r].name, None)
            del self.datasets[r]
        self._refresh()

    def _refresh(self):
        self.table.setRowCount(len(self.datasets))
        for i, d in enumerate(self.datasets):
            dts = np.unique(np.concatenate([f.dT for f in d.files]).round(2)) if d.files else []
            self.table.setItem(i, 0, QTableWidgetItem(d.name))
            self.table.setItem(i, 1, QTableWidgetItem(str(len(d.files))))
            txt = ", ".join(f"{x:g}" for x in dts[:6]) + (" ..." if len(dts) > 6 else "")
            self.table.setItem(i, 2, QTableWidgetItem(txt))
        cur = self.show_combo.currentText()
        self.show_combo.blockSignals(True)
        self.show_combo.clear()
        self.show_combo.addItems([d.name for d in self.datasets])
        if cur:
            self.show_combo.setCurrentText(cur)
        self.show_combo.blockSignals(False)

    # analysis
    def run(self):
        s = self.form_value("thermo")
        if s is None:
            return
        if not self.datasets:
            error_box(self, self.title, "Add a data set first.")
            return
        sets = list(self.datasets)

        def job():
            out, errs = {}, []
            for d in sets:
                try:
                    out[d.name] = th.analyse(d, s)
                except Exception as exc:  # noqa: BLE001
                    errs.append(f"{d.name}: {exc}")
            return out, errs

        def done(res):
            self.results, errs = res
            if errs:
                error_box(self, self.title, "\n".join(errs))
            lines = []
            for d in self.datasets:
                if d.name in self.results:
                    lines += self.results[d.name].summary_lines()
            self.result_label.setText("<br>".join(lines))
            self._draw()

        self.background(job, on_done=done, busy="Analysing thermovoltages ...")

    def current(self):
        return self.results.get(self.show_combo.currentText())

    def panels(self):
        r = self.current()
        return figs.thermo(r) if r is not None else []

    def _draw(self):
        self.plot.show_panels(self.panels(), 2)

    def has_results(self):
        return bool(self.results)

    # export
    def figure_prefix(self, prefix):
        r = self.current()
        return f"{prefix}_Thermo_{_safe(r.name)}" if r is not None else prefix

    def export_figures(self, folder, prefix):
        files = []
        for name, r in self.results.items():
            files += style.export_panels(figs.thermo(r), os.path.join(folder, "figures"),
                                         f"{prefix}_Thermo_{_safe(name)}", self.main.figure_settings)
        return files

    def export(self, folder, prefix):
        files = []
        for name, r in self.results.items():
            tag = f"{prefix}_Thermo_{_safe(name)}"
            f1 = os.path.join(folder, f"{tag}_peaks.txt")
            ex.save_columns(f1, ["dT_K", "V_peak_uV", "V_sigma_uV", "V_peak_err_uV", r.settings.unit],
                            r.dT_levels, r.v_peak, r.v_sigma, r.v_err, r.counts)
            f2 = os.path.join(folder, f"{tag}_seebeck_values.txt")
            cols = [r.S_values] + ([r.S_logG] if r.S_logG is not None else [])
            ex.save_columns(f2, ["S_uV_per_K"] + (["logG_G0"] if r.S_logG is not None else []), *cols)
            f3 = os.path.join(folder, f"{tag}_voltage_histograms.txt")
            names, out = [], []
            for d in r.dT_levels:
                _, _, _, h = th.gauss_peak(r.groups[d], r.settings.v_bins)
                if h is not None:
                    names += [f"V_uV_dT{d:g}", f"counts_dT{d:g}"]
                    out += [h[0], h[1]]
            if out:
                ex.save_columns(f3, names, *out)
                files.append(f3)
            f4 = ex.save_summary(folder, [l.replace("<b>", "").replace("</b>", "") for l in r.summary_lines()]
                                 + ["settings: " + str(r.settings.to_dict())], f"{tag}_summary.txt")
            files += [f1, f2, f4]
        return files


def _safe(name: str) -> str:
    return "".join(c if c.isalnum() or c in "-_." else "_" for c in name)
