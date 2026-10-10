"""Ensemble (EGaIn / large-area) junctions tab."""

from __future__ import annotations

import os

import numpy as np
from PySide6.QtWidgets import (
    QComboBox, QFileDialog, QFormLayout, QGroupBox, QInputDialog, QLabel, QTableWidget, QTableWidgetItem,
)

from .. import bj_io, egain as eg, export as ex, figures as figs, ivmodels as ivm, style
from .base import AnalysisTab
from .widgets import PlotPanel, SettingsForm, button_row, error_box, make_button, scroll_panel

EGAIN_LABELS = {"v_column": "bias column", "current_column": "current column", "current_unit": "current unit",
                "contact_diameter_um": "contact diameter (um)", "contact_area_cm2": "contact area (cm2, overrides)",
                "short_threshold": "short if |J| reaches (A/cm2)", "v_step": "bias grid step (V)",
                "bins_per_decade": "log|J| bins / decade", "alpha": "CI: 1 - alpha with alpha",
                "v_report": "report at bias (V)", "fit_slm": "single-level fit of <log|J|>",
                "temperature": "temperature (K)"}
EGAIN_CHOICES = {"current_unit": list(eg.UNIT_SCALE)}
EGAIN_TIPS = {
    "current_column": "Column number (0 = first) or header name. Every file is one junction.",
    "current_unit": "A/cm2: the column already holds the current density J.",
    "contact_diameter_um": "J = I / (pi d^2 / 4). Use the measured geometric contact diameter.",
    "short_threshold": "A junction whose |J| ever reaches this is a short: it is counted against the yield and "
                       "left out of the statistics.",
    "alpha": "Confidence interval t(1 - alpha/2, n_j - 1) sigma / sqrt(n_j - 1) with n_j the number of junctions "
             "(Reus et al. 2012).",
    "v_report": "Bias for the log|J| and rectification histograms and for beta.",
}


class EGaInTab(AnalysisTab):
    title = "EGaIn (ensemble)"

    def __init__(self, project, main):
        super().__init__(project, main)
        self.form = SettingsForm("Large-area junctions", eg.EGaInSettings(), EGAIN_LABELS, EGAIN_CHOICES,
                                 EGAIN_TIPS)
        self.forms["egain"] = self.form
        self.datasets: list[eg.Dataset] = []
        self.results: dict = {}
        self.beta = None
        self.table = QTableWidget(0, 3)
        self.table.setHorizontalHeaderLabels(["data set", "junctions", "length"])
        self.table.setMinimumHeight(140)
        self.table.itemChanged.connect(self._length_edited)
        self.show_combo = QComboBox()
        self.show_combo.currentIndexChanged.connect(lambda _: self._draw())
        box = QGroupBox("Data sets (one file per junction)")
        lay = QFormLayout(box)
        lay.addRow(button_row(make_button("Add files ...", self.add_files),
                              make_button("Add folder ...", self.add_folder)))
        lay.addRow(button_row(make_button("Add simulated", self.add_simulated,
                                          "Fake junctions with known <log|J|>, sigma, beta and rectification."),
                              make_button("Remove", self.remove)))
        lay.addRow(self.table)
        lay.addRow("show", self.show_combo)
        self.result_label = QLabel("Load J(V) or I(V) files of single junctions (EGaIn, CP-AFM, ...). "
                                   "Enter a molecular length per data set to get beta.")
        self.result_label.setWordWrap(True)
        self.layout_.addWidget(scroll_panel(box, self.form, make_button("Analyse all", self.run),
                                            make_button("Fit beta (log|J| vs length)", self.fit_beta),
                                            self.result_label, width=380))
        self.plot = PlotPanel()
        self.layout_.addWidget(self.plot, 1)

    # data sets
    def _add(self, name, junctions):
        if not junctions:
            return
        self.datasets.append(eg.Dataset(name, junctions))
        self._refresh()

    def _load(self, paths, name):
        s = self.form_value("egain")
        if s is None or not paths:
            return
        junctions, bad = [], []
        for p in paths:
            try:
                junctions.append(eg.load_junction(p, s))
            except Exception as exc:  # noqa: BLE001
                bad.append(f"{os.path.basename(p)}: {exc}")
        if bad:
            error_box(self, self.title, "Some files could not be read:\n" + "\n".join(bad[:20]))
        self._add(name, junctions)

    def add_files(self):
        paths, _ = QFileDialog.getOpenFileNames(self, "Junction files (one per junction)", self.main.last_dir,
                                                "Data (*.txt *.csv *.dat *.tsv *.pxp *.ibw);;All files (*)")
        if not paths:
            return
        self.main.last_dir = os.path.dirname(paths[0])
        name, ok = QInputDialog.getText(self, "Data set", "Name:", text=os.path.basename(os.path.dirname(paths[0])))
        if ok:
            self._load(paths, name.strip() or f"set {len(self.datasets) + 1}")

    def add_folder(self):
        folder = QFileDialog.getExistingDirectory(self, "Folder with junction files", self.main.last_dir)
        if not folder:
            return
        self.main.last_dir = folder
        paths = [p for p in bj_io.list_data_files(folder) if not p.lower().endswith(".tdms")]
        self._load(paths, os.path.basename(folder))

    def add_simulated(self):
        s = self.form_value("egain")
        if s is None:
            return
        k = len(self.datasets)
        L = 10 + 2 * k
        self.datasets.append(eg.Dataset(f"sim C{L}", eg.simulate(length=L, seed=k, area_cm2=s.area()), float(L)))
        self._refresh()

    def remove(self):
        rows = sorted({i.row() for i in self.table.selectedIndexes()}, reverse=True)
        for r in rows:
            name = self.datasets[r].name
            del self.datasets[r]
            self.results.pop(name, None)
        self._refresh()

    def _refresh(self):
        self.table.blockSignals(True)
        self.table.setRowCount(len(self.datasets))
        for i, d in enumerate(self.datasets):
            self.table.setItem(i, 0, QTableWidgetItem(d.name))
            self.table.setItem(i, 1, QTableWidgetItem(str(len(d.junctions))))
            self.table.setItem(i, 2, QTableWidgetItem("" if d.length is None else f"{d.length:g}"))
        self.table.blockSignals(False)
        cur = self.show_combo.currentText()
        self.show_combo.blockSignals(True)
        self.show_combo.clear()
        self.show_combo.addItems([d.name for d in self.datasets])
        if cur:
            self.show_combo.setCurrentText(cur)
        self.show_combo.blockSignals(False)

    def _length_edited(self, item):
        if item.column() == 0 and item.row() < len(self.datasets):
            old = self.datasets[item.row()].name
            new = item.text().strip() or old
            self.datasets[item.row()].name = new
            if old in self.results:
                self.results[new] = self.results.pop(old)
                self.results[new].name = new
            self._refresh()
            return
        if item.column() != 2 or item.row() >= len(self.datasets):
            return
        try:
            v = float(item.text()) if item.text().strip() else None
        except ValueError:
            v = None
        d = self.datasets[item.row()]
        d.length = v
        if d.name in self.results:
            self.results[d.name].length = v

    # analysis
    def run(self):
        s = self.form_value("egain")
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
                    out[d.name] = eg.analyse(d, s)
                except Exception as exc:  # noqa: BLE001
                    errs.append(f"{d.name}: {exc}")
            return out, errs

        def done(res):
            self.results, errs = res
            self.beta = None
            if errs:
                error_box(self, self.title, "\n".join(errs))
            self._report()
            self._draw()

        self.background(job, on_done=done, busy="Analysing junctions ...")

    def fit_beta(self):
        s = self.form_value("egain")
        if s is None:
            return
        res = [self.results[d.name] for d in self.datasets if d.name in self.results]
        self.beta = eg.fit_beta(res, s.v_report)
        if self.beta is None:
            error_box(self, self.title, "Analyse at least two data sets with a length entered in the table.")
            return
        self._report()
        self._draw()

    def _report(self):
        lines = []
        for d in self.datasets:
            r = self.results.get(d.name)
            if r is not None:
                lines += r.summary_lines()
        b = self.beta
        if b:
            lines.append(f"<b>beta = {b['beta']:.3f} +/- {b['beta_err']:.3f} per length unit, log J0 = "
                         f"{b['logJ0']:.2f} +/- {b['logJ0_err']:.2f}</b> (at {b['v']:+.2f} V, {len(b['points'])} sets)")
        self.result_label.setText("<br>".join(lines))

    def current(self):
        return self.results.get(self.show_combo.currentText())

    def panels(self):
        r = self.current()
        if r is None:
            return []
        return figs.egain(r, self.beta, ivm.slm_current)

    def _draw(self):
        self.plot.show_panels(self.panels(), 3)

    def has_results(self):
        return bool(self.results)

    # export
    def figure_prefix(self, prefix):
        r = self.current()
        return f"{prefix}_EGaIn_{_safe(r.name)}" if r is not None else prefix

    def export_figures(self, folder, prefix):
        files = []
        for name, r in self.results.items():
            files += style.export_panels(figs.egain(r, self.beta, ivm.slm_current), os.path.join(folder, "figures"),
                                         f"{prefix}_EGaIn_{_safe(name)}", self.main.figure_settings)
        return files

    def export(self, folder, prefix):
        files = []
        for name, r in self.results.items():
            tag = f"{prefix}_EGaIn_{_safe(name)}"
            f1 = os.path.join(folder, f"{tag}_logJ.txt")
            ex.save_columns(f1, ["V", "mean_logJ", "sigma_logJ", "CI_halfwidth", "sweeps"],
                            r.grid, r.mu, r.sigma, r.ci, r.counts)
            f2 = os.path.join(folder, f"{tag}_sweeps_logJ.txt")
            ex.save_columns(f2, ["V"] + [f"sweep{i}_junction{j}" for i, j in enumerate(r.sweep_junction)],
                            r.grid, *r.logJ)
            vs = sorted(r.logR)
            f3 = os.path.join(folder, f"{tag}_logR.txt")
            ex.save_columns(f3, ["V", "mean_logR", "sigma_logR", "cycles"], vs,
                            [eg.gauss_fit(r.logR[v], 20)[0] for v in vs], [eg.gauss_fit(r.logR[v], 20)[1] for v in vs],
                            [len(r.logR[v]) for v in vs])
            n = max(len(r.vt_plus), 1)
            pad = lambda a: np.r_[a, np.full(n - len(a), np.nan)]  # noqa: E731
            f4 = os.path.join(folder, f"{tag}_Vtrans.txt")
            ex.save_columns(f4, ["Vtrans_plus", "Vtrans_minus"], pad(r.vt_plus), pad(r.vt_minus))
            f5 = ex.save_summary(folder, r.summary_lines() + ["settings: " + str(r.settings.to_dict())],
                                 f"{tag}_summary.txt")
            files += [f1, f2, f3, f4, f5]
        if self.beta:
            b = self.beta
            files.append(ex.save_summary(folder, [
                f"beta = {b['beta']:.5g} +/- {b['beta_err']:.3g} per length unit",
                f"log10 J0 = {b['logJ0']:.5g} +/- {b['logJ0_err']:.3g} (J0 in A/cm2) at {b['v']:+.3f} V",
                "data set\tlength\tmean_logJ\tsigma_logJ"] +
                [f"{p[3]}\t{p[0]:g}\t{p[1]:.4f}\t{p[2]:.4f}" for p in b["points"]], f"{prefix}_EGaIn_beta.txt"))
        return files


def _safe(name: str) -> str:
    return "".join(c if c.isalnum() or c in "-_." else "_" for c in name)
