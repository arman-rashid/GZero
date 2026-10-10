"""Histogram, plateau length, junction statistics, correlation and clustering tabs."""

from __future__ import annotations

import os
from dataclasses import dataclass, asdict

import numpy as np
from PySide6.QtWidgets import QComboBox, QFormLayout, QGroupBox, QLabel, QTableWidget, QTableWidgetItem, QVBoxLayout

from .. import analysis as an, clustering as cl, figures as figs, junctions as jn
from .. import export as ex
from .base import AnalysisTab
from .widgets import PlotPanel, SettingsForm, error_box, make_button, scroll_panel

HIST_LABELS = {"g_min": "log G min", "g_max": "log G max", "bins_per_decade": "bins / decade",
               "z_min": "distance min (nm)", "z_max": "distance max (nm)", "z_bins": "distance bins",
               "normalize": "normalization"}
HIST_CHOICES = {"normalize": ["per trace", "counts", "density"]}


class ClusterSelector(QComboBox):
    """'All selected traces' or one cluster; kept in sync with the project."""

    def __init__(self, project):
        super().__init__()
        self.project = project
        project.selection_changed.connect(self.refresh)
        project.traces_changed.connect(self.refresh)
        self.refresh()

    def refresh(self):
        cur = self.currentIndex()
        self.blockSignals(True)
        self.clear()
        self.addItem("all accepted traces")
        for k in range(self.project.n_clusters):
            n = len(self.project.selected(k))
            self.addItem(f"cluster {k + 1} ({n} traces)")
        self.setCurrentIndex(cur if 0 <= cur < self.count() else 0)
        self.blockSignals(False)

    def cluster(self):
        return self.currentIndex() - 1 if self.currentIndex() > 0 else None

    def traces(self):
        return self.project.selected(self.cluster())

    def label(self):
        return "all" if self.cluster() is None else f"cluster{self.cluster() + 1}"


def _data_box(selector):
    box = QGroupBox("Data")
    lay = QFormLayout(box)
    lay.addRow("traces", selector)
    return box


# Histograms

@dataclass
class FitSettings:
    fit_min: float = -5.0
    fit_max: float = -1.5
    n_peaks: int = 1
    baseline: bool = True
    cmap: str = "viridis"
    color_max_percentile: float = 99.5
    overlay_clusters: bool = True


class HistogramTab(AnalysisTab):
    title = "Histograms"

    def __init__(self, project, main):
        super().__init__(project, main)
        self.hist_form = SettingsForm("Histogram", an.HistSettings(), HIST_LABELS, HIST_CHOICES)
        self.fit_form = SettingsForm("Peak fit / display", FitSettings(),
                                     {"fit_min": "fit from (log G)", "fit_max": "fit to (log G)",
                                      "n_peaks": "number of Gaussians", "baseline": "constant baseline",
                                      "cmap": "2D colormap", "color_max_percentile": "colour max (percentile)",
                                      "overlay_clusters": "overlay clusters in 1D"},
                                     {"cmap": ["viridis", "magma", "inferno", "jet", "turbo", "Greys", "hot"]})
        self.forms.update({"histogram": self.hist_form, "fit": self.fit_form})
        self.selector = ClusterSelector(project)
        self.result_label = QLabel("")
        self.result_label.setWordWrap(True)
        self.layout_.addWidget(scroll_panel(
            _data_box(self.selector), self.hist_form, self.fit_form,
            make_button("Update histograms", self.update),
            make_button("Fit peaks", self.fit), self.result_label))
        self.plot = PlotPanel()
        self.layout_.addWidget(self.plot, 1)
        self.res = {}
        self.last_peaks = None
        project.traces_changed.connect(self.update)
        self.selector.currentIndexChanged.connect(self.update)

    def update(self, *_):
        hs = self.form_value("histogram")
        fs = self.form_value("fit")
        if hs is None or fs is None:
            return
        traces = self.selector.traces()
        x, y = an.hist1d(traces, hs)
        zc, gc, H = an.hist2d(traces, hs)
        self.res = {"x": x, "y": y, "zc": zc, "gc": gc, "H": H, "n": len(traces), "hs": hs,
                    "label": self.selector.label(), "peaks": self.res.get("peaks") if self.res.get("n") == len(traces) else None}
        self._draw(fs)

    def fit(self):
        if not self.res:
            self.update()
        fs = self.form_value("fit")
        if fs is None or not self.res or self.res["n"] == 0:
            error_box(self, self.title, "No traces to fit.")
            return
        try:
            peaks, xf, yf, base = an.fit_peaks(self.res["x"], self.res["y"], fs.fit_min, fs.fit_max,
                                               fs.n_peaks, baseline=fs.baseline)
        except Exception as exc:  # noqa: BLE001
            error_box(self, "Peak fit failed", exc)
            return
        self.res["peaks"] = (peaks, xf, yf, base)
        self.last_peaks = (peaks, self.res["label"])   # kept when the view changes
        lines = []
        for k, p in enumerate(peaks):
            lines.append(f"Peak {k + 1}: log G = {p.center:.3f} +/- {p.center_err:.3f}  "
                         f"(G = {p.G:.3g} G0), sigma = {p.sigma:.3f}, FWHM = {p.fwhm:.3f} dec")
        self.result_label.setText("\n".join(lines))
        self._draw(fs)

    def _cluster_curves(self, hs, n_all):
        out = []
        for k in range(self.project.n_clusters):
            tk = self.project.selected(k)
            if tk:
                xk, yk = an.hist1d(tk, hs)
                scale = len(tk) / max(n_all, 1) if hs.normalize == "per trace" else 1.0
                out.append((f"Cluster {k + 1} ({len(tk)})", xk, yk * scale))
        return out

    def panels(self):
        r, fs = self.res, self.form_value("fit")
        if not r or fs is None:
            return []
        clusters = None
        if fs.overlay_clusters and self.selector.cluster() is None and self.project.n_clusters:
            clusters = self._cluster_curves(r["hs"], r["n"])
        label = "All accepted" if r["label"] == "all" else r["label"].replace("cluster", "Cluster ")
        return [figs.histogram_1d(r["x"], r["y"], r["hs"].normalize, label, r["n"], r.get("peaks"), clusters),
                figs.histogram_2d(r["zc"], r["gc"], r["H"], r["hs"].normalize, fs.cmap, fs.color_max_percentile)]

    def export_panels(self):
        """1D and 2D histograms of all accepted traces and of every cluster, one file each."""
        fs = self.form_value("fit")
        hs = self.res["hs"]
        out = []
        groups = [("all", self.project.selected())]
        groups += [(f"cluster{k + 1}", self.project.selected(k)) for k in range(self.project.n_clusters)]
        for label, traces in groups:
            if not traces:
                continue
            x, y = an.hist1d(traces, hs)
            zc, gc, H = an.hist2d(traces, hs)
            peaks = self.res.get("peaks") if label == self.res["label"] else None
            clusters = self._cluster_curves(hs, len(traces)) if label == "all" and fs.overlay_clusters else None
            name = "All accepted" if label == "all" else label.replace("cluster", "Cluster ")
            p1 = figs.histogram_1d(x, y, hs.normalize, name, len(traces), peaks, clusters or None)
            p2 = figs.histogram_2d(zc, gc, H, hs.normalize, fs.cmap, fs.color_max_percentile)
            if label != "all":
                p1.name, p2.name = f"{label}_{p1.name}", f"{label}_{p2.name}"
            out += [p1, p2]
        return out

    def _draw(self, fs=None):
        self.plot.show_panels(self.panels(), 2, [1, 1.3])

    def has_results(self):
        return bool(self.res) and self.res.get("n", 0) > 0

    def export(self, folder, prefix):
        """Histograms of all selected traces and of every cluster, plus the shown fit/figure."""
        if not self.has_results():
            return []
        r = self.res
        hs = r["hs"]
        groups = [("all", self.project.selected())]
        groups += [(f"cluster{k + 1}", self.project.selected(k)) for k in range(self.project.n_clusters)]
        files = []
        for label, traces in groups:
            if not traces:
                continue
            tag = prefix if label == "all" else f"{prefix}_{label}"
            x, y = an.hist1d(traces, hs)
            zc, gc, H = an.hist2d(traces, hs)
            f1, f2, f3 = (os.path.join(folder, f"{tag}_{s}.txt")
                          for s in ("logHist", "3D_Histogram", "3D_hist_scales"))
            ex.save_columns(f1, ["logG_G0", "counts"], x, y)
            ex.save_matrix(f2, H.T)  # rows = log G, columns = distance (as plotted)
            n = max(len(zc), len(gc))
            pad = lambda a: np.r_[a, np.full(n - len(a), np.nan)]  # noqa: E731
            ex.save_columns(f3, ["distance_nm", "logG_G0"], pad(zc), pad(gc))
            files += [f1, f2, f3]
        tag = f"{prefix}_{r['label']}" if r["label"] != "all" else prefix
        if r.get("peaks"):
            peaks = r["peaks"][0]
            p = os.path.join(folder, f"{tag}_peak_fit.txt")
            ex.save_columns(p, ["center_logG", "center_err", "G_G0", "sigma_dec", "fwhm_dec", "amplitude"],
                            [q.center for q in peaks], [q.center_err for q in peaks], [q.G for q in peaks],
                            [q.sigma for q in peaks], [q.fwhm for q in peaks], [q.amplitude for q in peaks])
            files.append(p)
        return files


# Plateau length

@dataclass
class PlateauSettings:
    window1_lo: float = -4.5
    window1_hi: float = -2.5
    window2_lo: float | None = None
    window2_hi: float | None = None
    window3_lo: float | None = None
    window3_hi: float | None = None
    method: str = "span"
    max_length: float = 2.0
    bins: int = 60
    min_length: float = 0.0
    beta_fit_lo: float = -5.0
    beta_fit_hi: float = -1.0
    beta_reference: float = 20.0


class PlateauTab(AnalysisTab):
    title = "Plateau length"

    def __init__(self, project, main):
        super().__init__(project, main)
        self.form = SettingsForm(
            "Plateau windows (log G)", PlateauSettings(),
            {"window1_lo": "window 1 from", "window1_hi": "window 1 to",
             "window2_lo": "window 2 from", "window2_hi": "window 2 to",
             "window3_lo": "window 3 from", "window3_hi": "window 3 to",
             "method": "length method", "max_length": "histogram max (nm)", "bins": "bins",
             "min_length": "ignore lengths below (nm)",
             "beta_fit_lo": "tunnelling fit from", "beta_fit_hi": "tunnelling fit to",
             "beta_reference": "expected beta (1/nm)"},
            {"method": ["span", "count"]},
            {"method": "span: from the first point below the upper limit to the last point above the lower limit.\n"
                       "count: number of points inside the window x distance step.",
             "beta_fit_lo": "log G range of the per-trace exponential tunnelling fit (used for the distance calibration).",
             "beta_reference": "Decay constant expected for your electrodes and medium (G ~ exp(-beta z)): "
                               "about 20 /nm for Au in vacuum (one decade per ~0.1 nm), lower in solvents."})
        self.forms["plateau"] = self.form
        self.selector = ClusterSelector(project)
        self.table = QTableWidget()
        self.table.setMinimumHeight(220)
        self.beta_label = QLabel("")
        self.beta_label.setWordWrap(True)
        self.layout_.addWidget(scroll_panel(_data_box(self.selector), self.form,
                                            make_button("Calculate", self.update), self.table,
                                            self.beta_label, width=360))
        self.plot = PlotPanel()
        self.layout_.addWidget(self.plot, 1)
        self.res = {}
        project.traces_changed.connect(self.update)
        self.selector.currentIndexChanged.connect(self.update)

    def _windows(self, ps):
        out = []
        for k in (1, 2, 3):
            lo, hi = getattr(ps, f"window{k}_lo"), getattr(ps, f"window{k}_hi")
            if lo is not None and hi is not None and hi > lo:
                out.append((lo, hi))
        return out

    def update(self, *_):
        ps = self.form_value("plateau")
        if ps is None:
            return
        traces = self.selector.traces()
        wins = self._windows(ps)
        lengths, stats = [], []
        for lo, hi in wins:
            L = an.plateau_lengths(traces, lo, hi, ps.method)
            lengths.append(L)
            stats.append(an.describe(L[L > ps.min_length]))
        beta = an.tunnelling_decay(traces, ps.beta_fit_lo, ps.beta_fit_hi) if traces else None
        self.res = {"windows": wins, "lengths": lengths, "stats": stats, "beta": beta, "ps": ps,
                    "label": self.selector.label(), "n": len(traces)}
        self._fill_table()
        self._draw()

    def _fill_table(self):
        keys = ["n", "mean", "median", "std", "skewness", "kurtosis", "geometric_mean", "min", "max"]
        r = self.res
        self.table.setRowCount(len(keys))
        self.table.setColumnCount(len(r["windows"]))
        self.table.setVerticalHeaderLabels(keys)
        self.table.setHorizontalHeaderLabels([f"[{lo:g}, {hi:g}]" for lo, hi in r["windows"]])
        for j, s in enumerate(r["stats"]):
            for i, k in enumerate(keys):
                v = s.get(k, np.nan)
                self.table.setItem(i, j, QTableWidgetItem(f"{v:.4g}" if isinstance(v, float) else str(v)))
        b = r["beta"]
        if b and b["n"]:
            ref = r["ps"].beta_reference
            self.beta_label.setText(
                f"Tunnelling decay: median beta = {b['beta_median']:.3g} per nm of the current distance axis "
                f"({b['n']} traces). To match the expected beta = {ref:g} /nm, multiply the displacement "
                f"ratio by {b['beta_median'] / ref:.3g}.")
        else:
            self.beta_label.setText("")

    def panels(self):
        r = self.res
        if not r:
            return []
        b = r["beta"]
        betas = b["betas"] if b and b["n"] else None
        return figs.plateau(r["windows"], r["lengths"], r["ps"].max_length, r["ps"].bins, r["ps"].min_length,
                            betas, b["beta_median"] if betas is not None else None)

    def _draw(self):
        self.plot.show_panels(self.panels(), 2, [2, 1])

    def figure_prefix(self, prefix):
        return f"{prefix}_{self.res['label']}" if self.res and self.res["label"] != "all" else prefix

    def has_results(self):
        return bool(self.res) and self.res["n"] > 0

    def export(self, folder, prefix):
        if not self.has_results():
            return []
        r = self.res
        tag = f"{prefix}_{r['label']}" if r["label"] != "all" else prefix
        f1 = os.path.join(folder, f"{tag}_plateau length.txt")
        ex.save_columns(f1, [f"length_nm_[{lo:g},{hi:g}]" for lo, hi in r["windows"]], *r["lengths"])
        f2 = os.path.join(folder, f"{tag}_Plateau_length_parameters.txt")
        lines = []
        for (lo, hi), s in zip(r["windows"], r["stats"]):
            lines.append(f"window [{lo:g}, {hi:g}] log G")
            lines += [f"  {k} = {v:.6g}" if isinstance(v, float) else f"  {k} = {v}" for k, v in s.items()]
        if r["beta"] and r["beta"]["n"]:
            lines.append(f"tunnelling beta median = {r['beta']['beta_median']:.6g} 1/nm ({r['beta']['n']} traces)")
        ex.save_summary(folder, lines, os.path.basename(f2))
        return [f1, f2]


# Junction statistics

JUNCTION_LABELS = {"g_lo": "molecular window from (log G)", "g_hi": "molecular window to (log G)",
                   "max_slope": "flat if slope below (dec/nm)", "slope_window": "slope over (nm)",
                   "min_plateau": "junction if flat length above (nm)", "block_size": "traces per block",
                   "g_bins": "conductance bins"}
JUNCTION_TIPS = {
    "max_slope": "Molecular plateaus fall well under 1 decade/nm, tunnelling about 4-9 decades/nm "
                 "(beta / ln 10). Points flatter than this count as plateau.",
    "slope_window": "Distance over which the local slope is measured; longer is smoother.",
    "min_plateau": "A trace counts as a molecular junction if it is flat inside the window for at least "
                   "this distance.",
    "block_size": "Consecutive traces grouped for the yield and conductance against trace number.",
}


class JunctionTab(AnalysisTab):
    title = "Junction statistics"

    def __init__(self, project, main):
        super().__init__(project, main)
        self.form = SettingsForm("Plateau detection", jn.JunctionSettings(), JUNCTION_LABELS, tips=JUNCTION_TIPS)
        self.forms["junctions"] = self.form
        self.selector = ClusterSelector(project)
        self.result_label = QLabel("")
        self.result_label.setWordWrap(True)
        info = QLabel("Per trace: is there a flat molecular plateau, how long is it and at what conductance. "
                      "The junction yield is the fraction of traces with a plateau (Kamenetska et al., "
                      "PRL 2009); against trace number it shows drifts during the measurement.")
        info.setWordWrap(True)
        self.layout_.addWidget(scroll_panel(_data_box(self.selector), self.form,
                                            make_button("Use peak from histogram fit", self.peak_from_histogram,
                                                        "Sets the window to the fitted peak centre +/- 3 sigma."),
                                            make_button("Calculate", self.update), self.result_label, info))
        self.plot = PlotPanel()
        self.layout_.addWidget(self.plot, 1)
        self.res = None
        self.label = "all"

    def peak_from_histogram(self):
        peaks = self.main.tab("Histograms").last_peaks
        if not peaks:
            error_box(self, self.title, "Fit the molecular peak in the Histograms tab first (Fit peaks).")
            return
        cur = self.form_value("junctions")
        if cur is None:
            return
        p = max(peaks[0], key=lambda q: q.amplitude)
        cur.g_lo, cur.g_hi = round(p.center - 3 * p.sigma, 3), round(p.center + 3 * p.sigma, 3)
        self.form.set_value(cur)
        self.status(f"Window set to [{cur.g_lo:g}, {cur.g_hi:g}]")

    def update(self):
        s = self.form_value("junctions")
        if s is None:
            return
        if s.g_hi <= s.g_lo:
            error_box(self, self.title, "The window must go from low to high log G.")
            return
        traces = self.selector.traces()
        if not traces:
            error_box(self, self.title, "No accepted traces. Detect the traces in the Data & Traces tab first.")
            return
        self.label = self.selector.label()
        self.background(jn.analyse, traces, s, on_done=self._done, busy="Looking for plateaus ...")

    def _done(self, res):
        self.res = res
        self.result_label.setText("<br>".join(res.summary_lines()))
        self.plot.show_panels(self.panels(), 3)

    def panels(self):
        return figs.junctions(self.res) if self.res is not None else []

    def figure_prefix(self, prefix):
        return f"{prefix}_{self.label}" if self.res is not None and self.label != "all" else prefix

    def has_results(self):
        return self.res is not None and self.res.n > 0

    def export(self, folder, prefix):
        if not self.has_results():
            return []
        r = self.res
        tag = self.figure_prefix(prefix)
        f1 = os.path.join(folder, f"{tag}_junctions.txt")
        ex.save_columns(f1, ["trace", "junction", "plateau_length_nm", "plateau_logG_G0"],
                        r.index, r.junction.astype(int), r.length, r.logG)
        bl = r.blocks()
        f2 = os.path.join(folder, f"{tag}_junction_blocks.txt")
        ex.save_columns(f2, ["mean_trace", "yield", "yield_se", "median_logG", "traces"], *bl.T)
        f3 = ex.save_summary(folder, r.summary_lines() + ["settings: " + str(r.settings.to_dict())],
                             f"{tag}_junction_summary.txt")
        return [f1, f2, f3]


# 2D cross-correlation

@dataclass
class CorrelationSettings:
    g_min: float = -6.0
    g_max: float = 0.5
    bins_per_decade: int = 20
    color_limit: float = 0.3
    cmap: str = "RdBu_r"


class CorrelationTab(AnalysisTab):
    title = "Correlation"

    def __init__(self, project, main):
        super().__init__(project, main)
        self.form = SettingsForm("2D cross-correlation", CorrelationSettings(),
                                 {"g_min": "log G min", "g_max": "log G max",
                                  "bins_per_decade": "bins / decade", "color_limit": "colour limit (+/-)",
                                  "cmap": "colormap"},
                                 {"cmap": ["RdBu_r", "seismic", "coolwarm", "PuOr_r", "jet"]})
        self.forms["correlation"] = self.form
        self.selector = ClusterSelector(project)
        info = QLabel("Pearson correlation of the per-trace histogram counts for every pair of "
                      "conductance bins (Makk et al., ACS Nano 2012). Positive off-diagonal "
                      "regions mark conductance values that occur together in the same trace.")
        info.setWordWrap(True)
        self.layout_.addWidget(scroll_panel(_data_box(self.selector), self.form,
                                            make_button("Calculate", self.update), info))
        self.plot = PlotPanel()
        self.layout_.addWidget(self.plot, 1)
        self.res = {}

    def update(self):
        s = self.form_value("correlation")
        if s is None:
            return
        traces = self.selector.traces()
        if len(traces) < 3:
            error_box(self, self.title, "At least 3 accepted traces are needed.")
            return
        n = max(4, int(round((s.g_max - s.g_min) * s.bins_per_decade)))
        edges = np.linspace(s.g_min, s.g_max, n + 1)
        C = an.correlation2d(traces, edges)
        h = an.hist_matrix(traces, edges).mean(axis=0)
        self.res = {"C": C, "edges": edges, "h": h, "label": self.selector.label(), "n": len(traces), "s": s}
        self.plot.show_panels(self.panels(), 2, [1.5, 1])

    def panels(self):
        r = self.res
        if not r:
            return []
        c = an.centers(r["edges"])
        return [figs.correlation_map(c, r["C"], r["s"].color_limit, r["s"].cmap), figs.mean_histogram(c, r["h"])]

    def figure_prefix(self, prefix):
        return f"{prefix}_{self.res['label']}" if self.res and self.res["label"] != "all" else prefix

    def has_results(self):
        return bool(self.res)

    def export(self, folder, prefix):
        if not self.res:
            return []
        r = self.res
        tag = f"{prefix}_{r['label']}" if r["label"] != "all" else prefix
        f1 = os.path.join(folder, f"{tag}_Correlation.txt")
        ex.save_matrix(f1, r["C"])
        f2 = os.path.join(folder, f"{tag}_Correlation_scales.txt")
        ex.save_columns(f2, ["logG_G0", "mean_counts_per_trace"], an.centers(r["edges"]), r["h"])
        return [f1, f2]


# Clustering

CLUSTER_LABELS = {
    "feature": "feature", "g_min": "log G min", "g_max": "log G max", "g_bins": "log G bins",
    "z_min": "distance min (nm)", "z_max": "distance max (nm)", "z_bins": "distance bins",
    "resample_points": "resample points", "pca_components": "PCA components (0 = none)",
    "method": "algorithm", "n_clusters": "number of clusters", "random_state": "random seed",
    "sort_by": "order clusters by", "k_scan": "k range for scan",
}
CLUSTER_CHOICES = {"feature": ["2D histogram", "1D histogram", "resampled trace"],
                   "method": ["k-means", "gaussian mixture", "agglomerative", "spectral"],
                   "sort_by": ["conductance", "size"]}


class ClusterTab(AnalysisTab):
    title = "Clustering"

    def __init__(self, project, main):
        super().__init__(project, main)
        self.form = SettingsForm("PCA + clustering", cl.ClusterSettings(), CLUSTER_LABELS, CLUSTER_CHOICES,
                                 {"feature": "Each trace becomes a normalised 2D (distance x log G) histogram image, "
                                             "which needs no reference and does not depend on trace length.",
                                  "method": "k-means: fast, spherical clusters\n"
                                            "gaussian mixture: elliptical clusters\n"
                                            "agglomerative (Ward): hierarchical\n"
                                            "spectral: non-convex shapes (slow for >5000 traces)"})
        self.forms["clustering"] = self.form
        self.info = QLabel("The cluster labels are stored with the traces. Choose a cluster in the 'traces' box "
                           "of the Histograms, Plateau length and Correlation tabs.")
        self.info.setWordWrap(True)
        self.result_label = QLabel("")
        self.result_label.setWordWrap(True)
        self.layout_.addWidget(scroll_panel(self.form, make_button("Run clustering", lambda: self.run(False)),
                                            make_button("Run + scan k (silhouette)", lambda: self.run(True)),
                                            make_button("Clear clusters", self.clear),
                                            self.result_label, self.info))
        self.plot = PlotPanel()
        self.layout_.addWidget(self.plot, 1)
        self.res = None
        self.traces = []
        project.traces_changed.connect(self.clear)

    def clear(self):
        for t in self.project.traces:
            t.cluster = -1
        self.project.n_clusters = 0
        self.res = None
        self.result_label.setText("")
        self.project.selection_changed.emit()

    def run(self, scan):
        cs = self.form_value("clustering")
        if cs is None:
            return
        traces = self.project.selected()
        if len(traces) < max(3, cs.n_clusters):
            error_box(self, self.title, "Not enough accepted traces.")
            return
        self.background(cl.run, traces, cs, scan, on_done=lambda r: self._done(r, traces, cs),
                        busy="Clustering ...")

    def _done(self, res, traces, cs):
        for t in self.project.traces:
            t.cluster = -1
        for t, lab in zip(traces, res.labels):
            t.cluster = int(lab)
        self.project.n_clusters = res.n_clusters
        self.res, self.traces, self.cs = res, traces, cs
        sizes = np.bincount(res.labels, minlength=res.n_clusters)
        txt = [f"{len(traces)} traces in {res.n_clusters} clusters; silhouette = {res.silhouette:.3f}"]
        txt += [f"Cluster {k + 1}: {n} traces ({100 * n / len(traces):.1f} %)" for k, n in enumerate(sizes)]
        if res.k_scan:
            best = max(res.k_scan, key=res.k_scan.get)
            txt.append(f"Silhouette scan: highest at k = {best} (also consider the physical meaning).")
        self.result_label.setText("\n".join(txt))
        self.project.selection_changed.emit()
        self._draw()

    def panels(self):
        if self.res is None:
            return []
        res, traces, cs = self.res, self.traces, self.cs
        hs = an.HistSettings(g_min=cs.g_min, g_max=cs.g_max)
        hs2 = an.HistSettings(g_min=cs.g_min, g_max=cs.g_max, bins_per_decade=30,
                              z_min=cs.z_min, z_max=cs.z_max, z_bins=80)
        h1, h2 = [], []
        for j in range(res.n_clusters):
            tj = [traces[i] for i in np.flatnonzero(res.labels == j)]
            x, y = an.hist1d(tj, hs)
            zc, gc, H = an.hist2d(tj, hs2)
            h1.append((x, y, len(tj)))
            h2.append((zc, gc, H, len(tj)))
        return figs.clustering(res, h1, h2)

    def _draw(self):
        self.plot.show_panels(self.panels(), 4 if self.res.k_scan else 3)

    def has_results(self):
        return self.res is not None

    def export(self, folder, prefix):
        if self.res is None:
            return []
        res, traces = self.res, self.traces
        sizes = np.bincount(res.labels, minlength=res.n_clusters)
        lines = [f"Total number of traces       = {len(traces)}"]
        lines += [f"Number of traces in Cluster {k + 1}= {n}" for k, n in enumerate(sizes)]
        lines.append(f"silhouette = {res.silhouette:.4f}")
        lines.append("settings: " + str(asdict(self.cs)))
        files = [ex.save_summary(folder, lines, f"{prefix}_Clustering_Results.txt")]
        f = os.path.join(folder, f"{prefix}_cluster_labels.txt")
        ex.save_columns(f, ["trace", "cluster", "PC1", "PC2", "PC3"],
                        [t.index for t in traces], res.labels + 1,
                        *[res.scores[:, i] if res.scores.shape[1] > i else np.zeros(len(traces)) for i in range(3)])
        files.append(f)
        hs = an.HistSettings(g_min=self.cs.g_min, g_max=self.cs.g_max)
        cols, names = [], []
        for j in range(res.n_clusters):
            tj = [traces[i] for i in np.flatnonzero(res.labels == j)]
            x, y = an.hist1d(tj, hs)
            if not cols:
                cols.append(x)
                names.append("logG_G0")
            cols.append(y)
            names.append(f"cluster{j + 1}_counts_per_trace")
        f2 = os.path.join(folder, f"{prefix}_cluster_logHist.txt")
        ex.save_columns(f2, names, *cols)
        files.append(f2)
        return files
