"""Main window."""

from __future__ import annotations

import json
import os
import sys

if __name__ == "__main__" and not __package__:
    # Started as a plain script (e.g. PyCharm "Run" on this file): relative imports
    # need the package context, so re-launch as  python -m gzero.
    import runpy
    sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
    runpy.run_module("gzero", run_name="__main__")
    sys.exit()

from PySide6.QtCore import QSettings
from PySide6.QtGui import QAction, QKeySequence
from PySide6.QtWidgets import (
    QApplication, QFileDialog, QInputDialog, QLabel, QMainWindow, QMessageBox, QTabWidget,
)

from .. import export as ex
from .project import Project
from .tabs_data import ConvertTab, DataTab, SimulationTab
from .tabs_egain import EGaInTab
from .tabs_noise import EventsTab, IVTab, ModulationTab, NoiseTab
from .tabs_stats import ClusterTab, CorrelationTab, HistogramTab, JunctionTab, PlateauTab
from .widgets import error_box

APP_NAME = "GZero"
APP_VERSION = "1.4"

ABOUT = f"""<h3>{APP_NAME} {APP_VERSION}</h3>
<p>Conversion and statistical analysis of single-molecule break-junction measurements
(MCBJ and STM-BJ).</p>
<ul>
<li><b>Convert TDMS</b>: channel calibration (linear / diode polynomial), Ia/Ib combination,
series-resistor compensation, filters, .cvr export.</li>
<li><b>Data &amp; Traces</b>: .cvr/.tdms/.txt loading, piezo-based segmentation with hold detection,
cutting limits, zero set, trace browser with manual good/bad marking.</li>
<li><b>Histograms</b>: 1D log-G and 2D conductance-distance histograms, multi-Gaussian fits.</li>
<li><b>Plateau length</b>: up to 3 windows, statistics, tunnelling-decay calibration.</li>
<li><b>Junction statistics</b>: per-trace plateau detection, junction yield (Kamenetska 2009),
plateau conductance and its drift over the measurement.</li>
<li><b>Correlation</b>: 2D cross-correlation histograms (Makk et al. 2012).</li>
<li><b>Clustering</b>: PCA + k-means / Gaussian mixture / agglomerative / spectral on per-trace 2D
histograms (Lemmer 2016, Cabosart 2019).</li>
<li><b>Flicker noise</b>: noise power vs G with scaling exponent n (NP ~ G<sup>n</sup>;
n ~ 1 through-bond, n ~ 2 through-space): bell-shaped 2D Gaussian analysis (Adak et al.,
Nano Lett. 2015) and the robust ADF + Theil-Sen method (Morris et al., J. Phys. Chem. C 2025).</li>
<li><b>I-V</b>: sweeps, dI/dV, single-level model fit (Zotti 2010) and transition voltage
spectroscopy (Beebe 2006, Baldea 2012).</li>
<li><b>Events</b>: flickering and mechanical events in G(t) and conductance-time histograms of holds
(Rashid et al., JACS 2025).</li>
<li><b>EGaIn (ensemble)</b>: large-area junctions, Gaussian log|J| statistics with junction-based
confidence intervals, yield, rectification, V<sub>trans</sub>, single-level fit and beta from a length series
(Reus 2012; XMe, Hong group; GaussFit, Chiechi group).</li>
<li><b>Piezo modulation</b> (lock-in beta), <b>Simulation</b> (synthetic ground truth).</li>
</ul>"""


class MainWindow(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle(f"{APP_NAME} {APP_VERSION}")
        self.settings = QSettings("GZero", "GZero")
        if not self.settings.allKeys():
            # first start after the rename: take over what the old version stored
            old = QSettings("BJAnalysis", "BJAnalysis")
            for key in old.allKeys():
                self.settings.setValue(key, old.value(key))
        self.last_dir = self.settings.value("last_dir", os.path.expanduser("~"))
        self.project = Project()
        self.figure_settings = self._load_figure_settings()
        self.tabs = QTabWidget()
        self.setCentralWidget(self.tabs)
        self._tabs = {}
        for cls in (ConvertTab, DataTab, HistogramTab, PlateauTab, JunctionTab, CorrelationTab, ClusterTab,
                    NoiseTab, EventsTab, IVTab, ModulationTab, EGaInTab, SimulationTab):
            tab = cls(self.project, self)
            self._tabs[tab.title] = tab
            self.tabs.addTab(tab, tab.title)
        self.tabs.setCurrentWidget(self._tabs["Data & Traces"])

        self.status_label = QLabel("Ready")
        self.statusBar().addWidget(self.status_label, 1)
        self.project.status.connect(self.status_label.setText)
        self._menus()
        self.resize(1500, 950)

    # tabs
    def tab(self, title):
        return self._tabs[title]

    def show_tab(self, title):
        self.tabs.setCurrentWidget(self._tabs[title])

    # menus
    def _menus(self):
        m = self.menuBar().addMenu("&File")
        for text, slot, key in (
            ("Open data files ...", lambda: self.tab("Data & Traces").load_files(), QKeySequence.Open),
            ("Open data folder ...", lambda: self.tab("Data & Traces").load_folder(), None),
            (None, None, None),
            ("Save results ...", self.save_results, QKeySequence.Save),
            ("Send plots to Origin", self.send_to_origin, None),
            (None, None, None),
            ("Figure export settings ...", self.edit_figure_settings, None),
            (None, None, None),
            ("Save settings ...", self.save_settings, None),
            ("Load settings ...", self.load_settings, None),
            (None, None, None),
            ("Exit", self.close, QKeySequence.Quit),
        ):
            if text is None:
                m.addSeparator()
                continue
            a = QAction(text, self)
            a.triggered.connect(slot)
            if key is not None:
                a.setShortcut(key)
            m.addAction(a)
        h = self.menuBar().addMenu("&Help")
        a = QAction("User manual", self)
        a.setShortcut(QKeySequence.HelpContents)
        a.triggered.connect(self.open_manual)
        h.addAction(a)
        a = QAction("About / methods", self)
        a.triggered.connect(lambda: QMessageBox.about(self, "About", ABOUT))
        h.addAction(a)

    # settings
    def _load_figure_settings(self):
        from .. import style
        try:
            d = json.loads(self.settings.value("figure_settings", "{}"))
            return style.FigureSettings(**d)
        except (TypeError, ValueError):
            return style.FigureSettings()

    def edit_figure_settings(self):
        from PySide6.QtWidgets import QDialog, QDialogButtonBox, QVBoxLayout
        from .. import style
        from .widgets import SettingsForm
        dlg = QDialog(self)
        dlg.setWindowTitle("Figure export settings")
        form = SettingsForm(
            "Exported figures (Nature format, one plot per file)", self.figure_settings,
            {"plot_with": "plot with", "width": "figure width", "font": "font", "font_size": "font size (pt)",
             "line_width": "data line width (pt)", "axes_width": "axis line width (pt)",
             "tick_direction": "tick direction", "formats": "file formats", "dpi": "raster resolution (dpi)"},
            {"plot_with": style.PLOT_PROGRAMS, "width": list(style.WIDTHS_MM), "tick_direction": ["out", "in"]},
            {"plot_with": "matplotlib: figure files in the 'figures' folder.\n"
                          "Origin: the plotted data go into Origin workbooks and every plot becomes an Origin "
                          "graph of the same size and style; the Origin project is saved in the results folder.\n"
                          "Needs Windows, Origin 2021 or newer and the originpro package "
                          "(py -m pip install originpro).",
             "font_size": "Nature: 5-7 pt. Tick labels and legends are drawn 1 pt smaller.",
             "line_width": "Nature: at least 0.5 pt.",
             "formats": "Comma-separated list of pdf, png, tiff, svg and eps (PDF keeps the text editable).",
             "dpi": "Nature: at least 300 dpi; 600 dpi is recommended for line art."})
        lay = QVBoxLayout(dlg)
        lay.addWidget(form)
        buttons = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
        buttons.accepted.connect(dlg.accept)
        buttons.rejected.connect(dlg.reject)
        lay.addWidget(buttons)
        if dlg.exec() == QDialog.Accepted:
            try:
                self.figure_settings = form.value()
            except ValueError as exc:
                error_box(self, "Figure export settings", exc)
                return
            self.settings.setValue("figure_settings", json.dumps(self.figure_settings.to_dict()))

    def all_settings(self) -> dict:
        d = {title: tab.settings_dict() for title, tab in self._tabs.items()}
        d["Figures"] = self.figure_settings.to_dict()
        return d

    def save_settings(self):
        path, _ = QFileDialog.getSaveFileName(self, "Save settings", os.path.join(self.last_dir, "bj_settings.json"),
                                              "JSON (*.json)")
        if path:
            with open(path, "w") as fh:
                json.dump(ex._jsonable(self.all_settings()), fh, indent=2)
            self.project.status.emit(f"Settings saved to {path}")

    def load_settings(self):
        path, _ = QFileDialog.getOpenFileName(self, "Load settings", self.last_dir,
                                              "JSON (*.json);;All files (*)")
        if not path:
            return
        try:
            with open(path) as fh:
                d = json.load(fh)
        except (OSError, ValueError) as exc:
            error_box(self, "Load settings", exc)
            return
        # accept both a pure settings file and a results *_config.json
        d = d.get("settings", d)
        for title, tab in self._tabs.items():
            if title in d:
                tab.apply_settings(d[title])
        if isinstance(d.get("Figures"), dict):
            from .. import style
            try:
                self.figure_settings = style.FigureSettings(**d["Figures"])
            except TypeError:
                pass
        self.project.status.emit(f"Settings loaded from {path}")

    # results
    def save_results(self):
        if not self.project.traces and not any(t.has_results() for t in self._tabs.values()):
            error_box(self, "Save results", "Nothing to save yet.")
            return
        base = QFileDialog.getExistingDirectory(self, "Choose where to create the results folder", self.last_dir)
        if not base:
            return
        default = self._result_name()
        name, ok = QInputDialog.getText(self, "Results folder", "Folder name:", text=default)
        if not ok or not name.strip():
            return
        name = name.strip()
        overwrite = False
        if os.path.exists(os.path.join(base, name)):
            ans = QMessageBox.question(self, "Folder already exists",
                                       f"'{name}' already exists.\nOverwrite its contents? "
                                       "(No creates a new numbered folder.)",
                                       QMessageBox.Yes | QMessageBox.No | QMessageBox.Cancel)
            if ans == QMessageBox.Cancel:
                return
            overwrite = ans == QMessageBox.Yes
        fs = self.figure_settings
        if fs.origin_graphs() and not self._origin_ready():
            return
        folder = ex.prepare_folder(base, name, overwrite)
        files, errors, origin_groups = [], [], []
        for title, tab in self._tabs.items():
            if not tab.has_results():
                continue
            try:
                files += tab.export(folder, name)
                if fs.figure_files():
                    files += tab.export_figures(folder, name)
                if fs.origin_graphs():
                    origin_groups.append((title, tab.figure_prefix(name), tab.export_panels()))
            except Exception as exc:  # noqa: BLE001
                errors.append(f"{title}: {exc}")
        if origin_groups:
            project = os.path.join(folder, f"{name}.opju")
            try:
                self._send_origin(origin_groups, project)
                files.append(project)
            except Exception as exc:  # noqa: BLE001
                errors.append(f"Origin: {exc}")
        config = {
            "recordings": [{"name": r.name, "source": r.meta.get("source"), "fs": r.fs, "samples": r.n}
                           for r in self.project.recordings],
            "trace_summary": self.project.summary_text(),
            "settings": self.all_settings(),
        }
        files.append(ex.save_config(folder, name, config))
        msg = f"Saved {len(files)} files to {folder}"
        self.project.status.emit(msg)
        if errors:
            error_box(self, "Some results could not be saved", "\n".join(errors))
        else:
            QMessageBox.information(self, "Results saved", msg)

    # Origin
    def _origin_ready(self) -> bool:
        from .. import origin
        reason = origin.unavailable_reason()
        if reason:
            error_box(self, "Origin", reason)
        return reason is None

    def _send_origin(self, groups, project_path=None):
        """Rebuild the panels in Origin (COM: must run in the GUI thread); busy cursor + progress in the status bar."""
        from PySide6.QtCore import Qt
        from .. import origin

        def progress(done, total, graph):
            self.project.status.emit(f"Origin: {done}/{total} {graph}")
            QApplication.processEvents()

        QApplication.setOverrideCursor(Qt.WaitCursor)
        try:
            return origin.send([g for g in groups if g[2]], self.figure_settings, project_path, progress)
        finally:
            QApplication.restoreOverrideCursor()

    def send_to_origin(self):
        """Current plots of every tab with results into Origin, without saving any files."""
        groups = []
        for title, tab in self._tabs.items():
            if tab.has_results():
                try:
                    groups.append((title, tab.figure_prefix(self._result_name()), tab.export_panels()))
                except Exception as exc:  # noqa: BLE001
                    error_box(self, "Send plots to Origin", f"{title}: {exc}")
                    return
        if not any(g[2] for g in groups):
            error_box(self, "Send plots to Origin", "Nothing to plot yet.")
            return
        if not self._origin_ready():
            return
        try:
            names = self._send_origin(groups)
        except Exception as exc:  # noqa: BLE001
            error_box(self, "Send plots to Origin", exc)
            return
        self.project.status.emit(f"Sent {len(names)} plots to Origin")

    def _result_name(self) -> str:
        return (self.project.recordings[0].name if self.project.recordings else "BJ") + "_Analysis"

    def open_manual(self):
        from PySide6.QtCore import QUrl
        from PySide6.QtGui import QDesktopServices
        path = manual_path()
        if path is None:
            error_box(self, "User manual", "The manual (docs/User_Manual.pdf) was not found next to the program.")
            return
        QDesktopServices.openUrl(QUrl.fromLocalFile(path))

    def closeEvent(self, event):
        self.settings.setValue("last_dir", self.last_dir)
        super().closeEvent(event)


def resource_dir() -> str:
    """Folder with bundled resources: the PyInstaller bundle, else the project folder."""
    if getattr(sys, "frozen", False):
        return getattr(sys, "_MEIPASS", os.path.dirname(sys.executable))
    return os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))


def manual_path() -> str | None:
    for base in (resource_dir(), os.path.dirname(sys.executable)):
        p = os.path.join(base, "docs", "User_Manual.pdf")
        if os.path.exists(p):
            return p
    return None


def main(argv=None):
    from .. import style
    app = QApplication.instance() or QApplication(sys.argv if argv is None else argv)
    app.setApplicationName(APP_NAME)
    style.light_palette(app)
    style.apply_style()
    w = MainWindow()
    w.show()
    return app.exec()


if __name__ == "__main__":
    sys.exit(main())
