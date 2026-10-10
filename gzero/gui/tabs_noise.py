"""Flicker noise, events, I-V and piezo modulation tabs."""

from __future__ import annotations

import os

import numpy as np
from PySide6.QtWidgets import QComboBox, QFormLayout, QGroupBox, QHBoxLayout, QLabel, QSpinBox, QVBoxLayout, QWidget

from .. import events as evm, export as ex, figures as figs, iv, ivmodels as ivm, modulation as md, noise as nz, style
from .base import AnalysisTab
from .widgets import PlotPanel, SettingsForm, decimate_for_plot, error_box, make_button, scroll_panel

NOISE_LABELS = {
    "source": "segments", "cut_initial_ms": "trim start (ms)", "cut_final_ms": "trim end (ms)",
    "window_ms": "window (ms, 0 = whole hold)",
    "max_windows_per_segment": "max windows / segment", "f_lo": "integrate from (Hz)",
    "f_hi": "integrate to (Hz)", "psd_method": "PSD method", "psd_window": "PSD window",
    "welch_nperseg": "Welch segment length",
    "detrend": "detrend", "g_min": "accept log G from", "g_max": "accept log G to",
    "max_abs_kurtosis": "max |kurtosis| (empty = off)", "max_drift": "max drift (decades, empty = off)",
    "drift_points": "drift averaging points",
    "stationarity_test": "ADF stationarity test", "adf_alpha": "ADF significance (alpha)",
    "peak_center": "peak centre (log G, empty = off)", "peak_sigma": "peak sigma (decades)",
    "peak_nsigma": "ends within +/- n sigma", "end_fraction": "end fraction",
    "end_points": "end points (0 = use fraction)",
    "subtract_floor": "subtract instrument floor",
    "floor_below": "floor = windows below log G", "estimator": "estimator for n",
    "n_min": "n scan from", "n_max": "n scan to",
    "n_step": "n step", "x_bins": "bell bins (log G)", "y_bins": "bell bins (log NP/G^n)",
    "bootstrap": "bootstrap samples",
}
NOISE_CHOICES = {"source": ["holds", "plateau windows"], "psd_method": ["welch", "periodogram"],
                 "psd_window": ["hann", "boxcar"], "detrend": ["constant", "linear"],
                 "estimator": ["Theil-Sen", "OLS", "2D Gaussian"]}
NOISE_TIPS = {
    "source": "holds: conductance-time traces recorded while the piezo is stopped (recommended)\n"
              "plateau windows: windows inside the selected breaking traces",
    "cut_initial_ms": "Time skipped after the piezo stops (mechanical relaxation); 5 ms in Morris et al. 2025.",
    "cut_final_ms": "Time skipped at the end of the hold; 5 ms in Morris et al. 2025.",
    "window_ms": "0: one noise value per trimmed hold (Morris et al. 2025). Above 0: split the holds into windows of this length.",
    "detrend": "constant: subtract the trace mean before the spectrum (this avoids an artificial 1/f^2 spectrum).",
    "max_abs_kurtosis": "Rejects windows with switching or spikes (excess kurtosis above this value).",
    "max_drift": "Rejects windows whose averaged log G at the start and at the end differ by more than this.",
    "stationarity_test": "Augmented Dickey-Fuller test on log10|G| (with a constant; lag chosen by AIC). A window is kept if p <= alpha.",
    "peak_center": "Molecular peak from the 1D histogram fit. The first and last 'end fraction' of each window "
                   "must average within centre +/- n sigma. Use 'Peak from histogram fit' to fill these fields.",
    "subtract_floor": "Subtracts the mean PSD of the windows below the floor level (amplifier noise).",
    "estimator": "Theil-Sen: median of the pairwise slopes; robust to outliers (Morris et al. 2025).\n"
                 "OLS: least squares, equivalent to zero Pearson correlation.\n"
                 "2D Gaussian: zero correlation of the fitted bell (Adak et al. 2015).",
}


def _recording_box(tab, combo):
    box = QGroupBox("Recording")
    lay = QFormLayout(box)
    lay.addRow("file", combo)
    tab.project.recordings_changed.connect(lambda: _fill_recordings(tab.project, combo))
    _fill_recordings(tab.project, combo)
    return box


def _fill_recordings(project, combo):
    cur = combo.currentText()
    combo.clear()
    combo.addItems([r.name for r in project.recordings])
    if cur:
        combo.setCurrentText(cur)


# Flicker noise

class NoiseTab(AnalysisTab):
    title = "Flicker noise"

    def __init__(self, project, main):
        super().__init__(project, main)
        self.form = SettingsForm("Noise power and selection", nz.NoiseSettings(),
                                 NOISE_LABELS, NOISE_CHOICES, NOISE_TIPS)
        self.forms["noise"] = self.form
        self.result_label = QLabel("Detect traces first; holds are found automatically "
                                   "when the piezo stops.")
        self.result_label.setWordWrap(True)
        self.seg_idx = QSpinBox()
        self.seg_idx.setRange(0, 0)
        self.seg_idx.valueChanged.connect(lambda _: self._draw())
        seg_row = QWidget()
        sl = QHBoxLayout(seg_row)
        sl.setContentsMargins(0, 0, 0, 0)
        sl.addWidget(QLabel("show segment"))
        sl.addWidget(self.seg_idx)
        presets = QGroupBox("Method presets")
        pl = QVBoxLayout(presets)
        pl.addWidget(make_button("Robust: Morris et al. 2025 (ADF + Theil-Sen)", self.preset_morris,
                                 "Trims 5 ms at both ends, uses the whole hold and a mean-subtracted one-sided spectrum "
                                 "(100-1000 Hz), applies the ADF stationarity test (alpha 0.05), keeps ends within "
                                 "+/- 2 sigma of the peak, and uses the Theil-Sen estimator."))
        pl.addWidget(make_button("Classic: Adak et al. 2015 (2D Gaussian bell)", self.preset_adak,
                                 "Uses 100 ms windows and a Welch PSD (100-1000 Hz); n is where the fitted 2D Gaussian "
                                 "has zero correlation."))
        pl.addWidget(make_button("Rashid et al. 2025 (Pearson r scan)", self.preset_rashid,
                                 "Drops the first 10 ms of each hold and uses the rest as one window, keeps holds whose "
                                 "first and last 100 points average within +/- 1 sigma of the peak, DFT squared "
                                 "(periodogram) integrated 100-1000 Hz, n where the Pearson r is smallest "
                                 "(scan 0.3-2.3 in steps of 0.01)."))
        pl.addWidget(make_button("Peak from histogram fit", self.peak_from_histogram,
                                 "Copies the centre and sigma of the fitted 1D-histogram peak that lies inside the accepted "
                                 "conductance range (Histograms tab, Fit peaks)."))
        info = QLabel(
            "NP ~ G<sup>n</sup>: <b>n ~ 1 through-bond</b>, <b>n ~ 2 through-space</b> coupling. "
            "The bell plot is the 2D histogram of log(NP/G<sup>n</sup>) vs log G at the chosen n; "
            "for the right n it is uncorrelated. Always report n with its standard error.")
        info.setWordWrap(True)
        self.layout_.addWidget(scroll_panel(presets, self.form, make_button("Run noise analysis", self.run),
                                            seg_row, self.result_label, info, width=380))
        self.plot = PlotPanel()
        self.layout_.addWidget(self.plot, 1)
        self.res = None
        self.scaling = None
        self.segments = []

    # presets
    def preset_morris(self):
        cur = self.form_value("noise")
        if cur is not None:
            self.form.set_value(nz.morris2025_settings(cur))
            self.status("Morris et al. 2025 settings applied. Set the peak (Peak from histogram fit), then run.")

    def preset_rashid(self):
        cur = self.form_value("noise")
        if cur is not None:
            self.form.set_value(nz.rashid2025_settings(cur))
            self.status("Rashid et al. 2025 settings applied. Set the peak (Peak from histogram fit), then run.")

    def preset_adak(self):
        cur = self.form_value("noise")
        if cur is None:
            return
        s = nz.NoiseSettings(g_min=cur.g_min, g_max=cur.g_max, source=cur.source)
        s.estimator = "2D Gaussian"
        self.form.set_value(s)
        self.status("Adak et al. 2015 settings applied.")

    def peak_from_histogram(self):
        hist = self.main.tab("Histograms")
        peaks = hist.last_peaks
        if not peaks:
            error_box(self, self.title, "Fit the molecular peak in the Histograms tab first (Fit peaks).")
            return
        cur = self.form_value("noise")
        if cur is None:
            return
        cands = [p for p in peaks[0] if cur.g_min <= p.center <= cur.g_max] or list(peaks[0])
        p = max(cands, key=lambda q: q.amplitude)
        cur.peak_center, cur.peak_sigma = round(p.center, 4), round(p.sigma, 4)
        self.form.set_value(cur)
        self.status(f"Peak set: log G = {p.center:.3f}, sigma = {p.sigma:.3f}")

    # run
    def _segments(self, ns):
        if ns.source == "holds":
            return self.project.holds
        return self.project.selected()

    def run(self):
        ns = self.form_value("noise")
        if ns is None:
            return
        segs = self._segments(ns)
        if not segs:
            error_box(self, self.title, "No hold segments were found. Detect the traces in the Data & Traces tab first "
                                        "(holds need a piezo signal), or set 'segments' to 'plateau windows'.")
            return

        def job():
            res = nz.analyse(segs, ns)
            G, NP = res.arrays()
            sc = None
            err = None
            if len(G) >= 5:
                try:
                    sc = nz.scaling_exponent(G, NP, ns)
                except Exception as exc:  # noqa: BLE001
                    err = str(exc)
            else:
                err = (f"Only {len(G)} accepted window{'' if len(G) == 1 else 's'}. "
                       "Relax the selection criteria.")
            return res, sc, err

        def done(out):
            self.res, self.scaling, err = out
            self.segments = segs
            self.ns = ns
            self.seg_idx.blockSignals(True)
            self.seg_idx.setRange(0, max(0, len(segs) - 1))
            self.seg_idx.blockSignals(False)
            self._report(err)
            self._draw()

        self.background(job, on_done=done, busy="Computing noise power ...")

    def _report(self, err):
        r = self.res
        acc = r.accepted()
        lines = [f"{len(r.windows)} windows from {len(self.segments)} segments, {len(acc)} accepted"]
        rej = r.rejected_counts()
        if rej:
            lines.append("Rejected: " + ", ".join(f"{k}: {v}" for k, v in rej.items()))
        lines.append(f"Mean PSD ~ 1/f^alpha with alpha = {r.alpha:.2f}")
        sc = self.scaling
        if sc:
            lo, hi = sc["n_tse_ci95"]
            lines += [f"<b>n (Theil-Sen) = {sc['n_tse']:.3f} &plusmn; {sc['n_tse_se']:.3f}</b> "
                      f"(95 % CI {lo:.3f} - {hi:.3f})",
                      f"n (OLS) = {sc['n_ols']:.3f} &plusmn; {sc['n_ols_se']:.3f}; "
                      f"minimum |Pearson r| on the scan grid at n = {sc['n_pearson_grid']:.2f}",
                      f"n (2D Gaussian fit) = {sc['n_fit']:.3f}",
                      f"Bell at n = {sc['n_best']:.3f} ({sc['estimator']}); centre at log G = {sc['bell'].x0:.2f}",
                      f"<b>{sc['interpretation']}</b>"]
        if err:
            lines.append(f"<span style='color:#b00'>{err}</span>")
        self.result_label.setText("<br>".join(lines))

    def panels(self):
        r = self.res
        if r is None:
            return []
        si = self.seg_idx.value()
        seg = self.segments[si] if 0 <= si < len(self.segments) else None
        panels = figs.noise(r, self.scaling, self.ns, seg, si)
        if not self.scaling:
            n_acc = len(r.accepted())

            def message(fig, ax):
                ax.text(0.5, 0.5, f"{n_acc} accepted window{'' if n_acc == 1 else 's'}:\n"
                        "not enough for a fit\n(relax the selection criteria)", ha="center", va="center",
                        transform=ax.transAxes, fontsize=9, color="0.4")
                ax.set_axis_off()
            panels.append(style.Panel("_message", message))
        return panels

    def export_panels(self):
        return [p for p in self.panels() if not p.name.startswith("_")]

    def _draw(self):
        self.plot.show_panels(self.panels(), 3)

    def has_results(self):
        return self.res is not None

    def export(self, folder, prefix):
        if self.res is None:
            return []
        r = self.res
        files = []
        f = os.path.join(folder, f"{prefix}_Noise_windows.txt")
        ws = r.windows
        ex.save_columns(f, ["segment", "start_sample", "samples", "G_G0", "noise_power_G0^2", "kurtosis",
                            "drift_dec", "adf_p", "accepted"],
                        [w.segment for w in ws], [w.start for w in ws], [w.length for w in ws],
                        [w.G for w in ws], [w.NP for w in ws], [w.kurtosis for w in ws],
                        [w.drift for w in ws], [w.adf_p for w in ws], [int(w.accepted) for w in ws])
        files.append(f)
        f = os.path.join(folder, f"{prefix}_Noise_PSD.txt")
        floor = r.psd_floor if r.psd_floor is not None else np.full_like(r.freqs, np.nan)
        ex.save_columns(f, ["frequency_Hz", "mean_PSD_G0^2/Hz", "floor_PSD"], r.freqs, r.psd_mean, floor)
        files.append(f)
        sc = self.scaling
        lines = [f"windows = {len(ws)}, accepted = {len(r.accepted())}",
                 f"rejected = {r.rejected_counts()}", f"alpha (1/f^alpha) = {r.alpha:.4f}"]
        if sc:
            f = os.path.join(folder, f"{prefix}_Noise_scaling.txt")
            ex.save_columns(f, ["n", "pearson_r", "kendall_tau", "gauss2d_rho"], sc["n_grid"],
                            sc["rho_pearson"], sc["tau_kendall"], sc["rho_fit"])
            files.append(f)
            f = os.path.join(folder, f"{prefix}_Noise_bell_hist.txt")
            ex.save_matrix(f, sc["hist"].T)
            files.append(f)
            b = sc["bell"]
            lines += [f"n_TheilSen = {sc['n_tse']:.4f} +/- {sc['n_tse_se']:.4f} (SE), "
                      f"95% CI {sc['n_tse_ci95'][0]:.4f} .. {sc['n_tse_ci95'][1]:.4f}, points {sc['tse_points']}",
                      f"n_OLS = {sc['n_ols']:.4f} +/- {sc['n_ols_se']:.4f} (SE), "
                      f"bootstrap 95% CI {sc['n_ols_ci95'][0]:.4f} .. {sc['n_ols_ci95'][1]:.4f}",
                      f"n_2DGaussian = {sc['n_fit']:.4f}",
                      f"n_Pearson_grid_minimum = {sc['n_pearson_grid']:.4f}",
                      f"estimator used for the bell = {sc['estimator']}, n = {sc['n_best']:.4f}",
                      f"bell: x0 = {b.x0:.4f}, y0 = {b.y0:.4f}, sx = {b.sx:.4f}, sy = {b.sy:.4f}, rho = {b.rho:.4f}",
                      f"bell x range (log G): {sc['x_edges'][0]:.4f} .. {sc['x_edges'][-1]:.4f}",
                      f"bell y range (log NP/G^n): {sc['y_edges'][0]:.4f} .. {sc['y_edges'][-1]:.4f}",
                      sc["interpretation"]]
        lines.append("settings: " + str(r.settings))
        files.append(ex.save_summary(folder, lines, f"{prefix}_Noise_summary.txt"))
        return files


# Flickering and mechanical events (Rashid et al. 2025, SI 2.6)

EVENT_LABELS = {"source": "segments", "lowpass_hz": "low-pass (Hz)", "lowpass_order": "low-pass order",
                "sg_side_points": "Savitzky-Golay side points", "sg_order": "Savitzky-Golay order",
                "dt_s": "derivative dt (s)", "min_width": "peak width (samples)",
                "peak_method": "peak detection", "polarity": "count",
                "flicker_threshold": "flicker threshold", "mechanical_threshold": "mechanical threshold",
                "g_min": "valid from log G", "g_max": "valid to log G", "level_ms": "level before/after (ms)",
                "t_bins": "hold histogram: time bins", "g_bins_per_decade": "hold histogram: bins / decade"}
EVENT_CHOICES = {"source": ["traces", "holds"], "peak_method": ["quadratic fit", "minimum width"],
                 "polarity": ["rises", "drops", "both"]}
EVENT_TIPS = {
    "dt_s": "The derivative is (G[i+1] - G[i-1]) / (2 dt) with G in G0. Rashid et al. 2025 used dt = 20 ms; the "
            "thresholds are on that scale.",
    "peak_method": "quadratic fit: like LabVIEW's Peak Detector, a quadratic is fitted over 'width' points and its "
                   "maximum must reach the threshold.\nminimum width: peaks at least 'width' samples wide at "
                   "half height.",
    "polarity": "rises: conductance increases only (LabVIEW Peak Detector default). both: also drops.",
    "g_min": "Events where the junction is outside this range (contact, noise floor) are false positives.",
}


class EventsTab(AnalysisTab):
    title = "Events"

    def __init__(self, project, main):
        super().__init__(project, main)
        self.form = SettingsForm("Flicker / mechanical event detection", evm.EventSettings(), EVENT_LABELS,
                                 EVENT_CHOICES, EVENT_TIPS)
        self.forms["events"] = self.form
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
        info = QLabel("G(t) is low-pass filtered (Butterworth) and smoothed (Savitzky-Golay), differentiated, and "
                      "peaks of the derivative above the thresholds count as flickers or mechanical events "
                      "(Rashid et al., JACS 2025). With holds, the conductance-time histogram of the holds is "
                      "shown too.")
        info.setWordWrap(True)
        self.layout_.addWidget(scroll_panel(self.form, make_button("Detect events", self.run), row,
                                            self.result_label, info, width=360))
        self.plot = PlotPanel()
        self.layout_.addWidget(self.plot, 1)
        self.res = None
        self.segs = []
        self.th = None

    def run(self):
        es = self.form_value("events")
        if es is None:
            return
        segs = self.project.holds if es.source == "holds" else self.project.selected()
        if not segs:
            error_box(self, self.title, "Nothing to analyse. Detect the traces in the Data & Traces tab first "
                                        "(holds need a piezo signal).")
            return

        def job():
            res = evm.analyse(segs, es)
            th = evm.time_histogram(segs, es) if es.source == "holds" else None
            return res, th

        def done(out):
            self.res, self.th = out
            self.segs = segs
            first = next((e.trace for e in self.res.events), 0)
            self.example.blockSignals(True)
            self.example.setRange(0, max(0, len(segs) - 1))
            self.example.setValue(first)
            self.example.blockSignals(False)
            self.result_label.setText("<br>".join(self.res.summary_lines()))
            self._draw()

        self.background(job, on_done=done, busy="Detecting events ...")

    def panels(self):
        if self.res is None:
            return []
        k = self.example.value()
        example = None
        if 0 <= k < len(self.segs):
            evs, yf, d, _ = evm.detect(self.segs[k], self.res.settings, k)
            kind = self.res.kinds[k] if k < len(self.res.kinds) else ""
            example = (self.segs[k], yf, d, evs, f"{self.res.settings.source[:-1]} {k}: {kind}")
        return figs.events(self.res, example, self.th)

    def _draw(self):
        self.plot.show_panels(self.panels(), 3)

    def has_results(self):
        return self.res is not None

    def export(self, folder, prefix):
        if self.res is None:
            return []
        r = self.res
        f1 = os.path.join(folder, f"{prefix}_events.txt")
        ex.save_columns(f1, ["segment", "sample", "time_s", "mechanical", "dGdt", "logG_before", "logG_after"],
                        [e.trace for e in r.events], [e.index for e in r.events],
                        [e.index / self.segs[e.trace].fs for e in r.events],
                        [int(e.kind == "mechanical") for e in r.events], [e.height for e in r.events],
                        [e.logG_before for e in r.events], [e.logG_after for e in r.events])
        f2 = os.path.join(folder, f"{prefix}_events_per_trace.txt")
        cls = {"quiet": 0, "flickering": 1, "mechanical": 2}
        ex.save_columns(f2, ["trace", "flickers", "mechanical", "class_0quiet_1flicker_2mechanical"],
                        r.index, r.n_flicker, r.n_mechanical, [cls[k] for k in r.kinds])
        files = [f1, f2]
        if self.th is not None:
            tc, gc, H, mp, sd = self.th
            f3 = os.path.join(folder, f"{prefix}_hold_time_histogram.txt")
            ex.save_matrix(f3, H.T)
            f4 = os.path.join(folder, f"{prefix}_hold_most_probable.txt")
            ex.save_columns(f4, ["time_s", "most_probable_logG", "sigma"], tc, mp, sd)
            n = max(len(tc), len(gc))
            pad = lambda a: np.r_[a, np.full(n - len(a), np.nan)]  # noqa: E731
            f5 = os.path.join(folder, f"{prefix}_hold_time_histogram_scales.txt")
            ex.save_columns(f5, ["time_s", "logG_G0"], pad(tc), pad(gc))
            files += [f3, f4, f5]
        files.append(ex.save_summary(folder, r.summary_lines() + ["settings: " + str(r.settings.to_dict())],
                                     f"{prefix}_events_summary.txt"))
        return files


# I-V

IV_LABELS = {"min_amplitude": "min sweep amplitude (V)", "smooth_points": "smoothing points",
             "poly_order": "polynomial order", "v_min": "curve must reach V below",
             "v_max": "curve must reach V above", "v_bins": "V bins", "i_bins": "I bins",
             "log_current": "log |I| histogram"}
IVM_LABELS = {"fit_v_max": "fit only |V| below (V, empty = all)", "fit_asymmetry": "fit asymmetry a",
              "n_molecules": "molecules in parallel N", "tvs_v_min": "TVS: ignore |V| below (V)",
              "max_curves": "fit at most N single curves", "temperature": "temperature (K, 0 = T = 0 form)",
              "integration": "mean-curve integration", "energy_step": "energy grid step (meV)",
              "energy_min": "energy grid from (eV)", "energy_max": "energy grid to (eV)",
              "fit_method": "fit algorithm", "fit_n": "fit N too (ensembles)",
              "mp_v_bins": "most probable: bias bins", "mp_i_bins": "most probable: current bins",
              "weight_by_sigma": "weight by 1/sigma"}
IVM_CHOICES = {"integration": ["analytic", "energy grid"], "fit_method": ["trust region", "Levenberg-Marquardt"]}
IVM_TIPS = {
    "fit_v_max": "The single-level model holds while the level is outside the bias window; "
                 "restrict the fit if the curves go much further.",
    "fit_asymmetry": "Lets the level move with the bias, e(V) = eps0 + a V. Off: symmetric junction (a = 0).",
    "n_molecules": "N in I = N G0 Gamma [...]. Keep 1 for single-molecule junctions.",
    "tvs_v_min": "Below this the current is too small for a stable V^2/|I|.",
    "temperature": "Fermi smearing of the leads. 0 uses the zero-temperature closed form; above 0 the exact "
                   "finite-temperature current (digamma form, or the energy grid below).",
    "integration": "analytic: closed form (fast). energy grid: Landauer integral on an energy grid as in "
                   "Rashid et al. 2025; used for the mean and most probable curves, the single curves always "
                   "use the closed form (both agree to ~1e-9).",
    "energy_step": "Rashid et al. 2025: 0.02 meV.",
    "fit_method": "trust region: bounded least squares. Levenberg-Marquardt: unbounded, as in Rashid et al. 2025.",
    "fit_n": "Lets the number of molecules N float, e.g. for EGaIn / large-area junctions. N and Gamma are "
             "strongly correlated; eps0 is the robust number then.",
    "mp_v_bins": "The most probable curve: in each bias bin a Gaussian is fitted to the current distribution "
                 "(Rashid et al. 2025, Fig. S22).",
}


class IVTab(AnalysisTab):
    title = "I-V"

    def __init__(self, project, main):
        super().__init__(project, main)
        self.rec = QComboBox()
        self.form = SettingsForm("I-V sweeps", iv.IVSettings(), IV_LABELS)
        self.forms["iv"] = self.form
        self.model_form = SettingsForm("Single-level model and TVS", ivm.IVModelSettings(), IVM_LABELS,
                                       IVM_CHOICES, IVM_TIPS)
        self.forms["iv_models"] = self.model_form
        self.result_label = QLabel("Requires a time-resolved bias: a TDMS file with a bias channel, or a text "
                                   "file with a bias column.")
        self.result_label.setWordWrap(True)
        self.model_label = QLabel("")
        self.model_label.setWordWrap(True)
        info = QLabel("Single-level model: one level at eps0 with coupling Gamma (Lorentzian transmission). "
                      "TVS: V<sub>t</sub> is the minimum of the Fowler-Nordheim plot; for this model "
                      "eps0 = (&radic;3/2) V<sub>t</sub> (Baldea 2012).")
        info.setWordWrap(True)
        self.layout_.addWidget(scroll_panel(_recording_box(self, self.rec), self.form,
                                            make_button("Analyse I-V", self.run), self.result_label,
                                            make_button("Preset: Rashid et al. 2025 (300 K, energy grid, LM)",
                                                        self.preset_rashid,
                                                        "Finite temperature (300 K), Landauer integral on an energy "
                                                        "grid from -10 to 10 eV in 0.02 meV steps, Levenberg-"
                                                        "Marquardt fit of the most probable I-V curve, "
                                                        "Gamma_L = Gamma_R."),
                                            self.model_form,
                                            make_button("Fit single-level model + TVS", self.run_models),
                                            self.model_label, info))
        self.plot = PlotPanel()
        self.layout_.addWidget(self.plot, 1)
        self.curves = []
        self.models = None

    def run(self):
        s = self.form_value("iv")
        rec = self.project.recording(self.rec.currentText())
        if s is None or rec is None:
            return
        self.background(iv.curves_from_recording, rec, s, on_done=lambda c: self._done(c, s),
                        busy="Finding I-V sweeps ...")

    def _done(self, curves, s):
        self.curves, self.s = curves, s
        self.models = None
        self.model_label.setText("")
        nf = sum(c.direction == "forward" for c in curves)
        if curves:
            self.result_label.setText(f"{len(curves)} curves ({nf} forward, {len(curves) - nf} backward)")
        else:
            self.result_label.setText("No bias sweeps were found in this recording (the bias is constant or "
                                      "varies less than the minimum sweep amplitude).")
        self._draw()

    def preset_rashid(self):
        cur = self.form_value("iv_models")
        if cur is None:
            return
        self.model_form.set_value(ivm.rashid2025_settings(cur))
        self.status("Rashid et al. 2025 settings applied (temperature 300 K assumed: the SI does not state it).")

    def run_models(self):
        ms = self.form_value("iv_models")
        if ms is None:
            return
        if not self.curves:
            error_box(self, self.title, "Find the I-V curves first (Analyse I-V).")
            return
        self.background(ivm.analyse, self.curves, ms, iv.mean_curve, on_done=self._models_done,
                        busy="Fitting the single-level model ...")

    def _models_done(self, res):
        self.models = res
        self.model_label.setText("<br>".join(res.summary_lines()) or "No fit converged.")
        self._draw()

    def _draw(self):
        self.plot.show_panels(self.panels(), 3 if self.models else 2)

    def panels(self):
        if not self.curves:
            return []
        panels = figs.iv(self.curves, self.s, iv.iv_hist2d, iv.mean_curve)
        if self.models is not None:
            panels += figs.iv_models(self.models, ivm.slm_current, ivm.fowler_nordheim)
        return panels

    def has_results(self):
        return bool(self.curves)

    def export(self, folder, prefix):
        if not self.curves:
            return []
        files = []
        sub = os.path.join(folder, f"{prefix}_IV_DATA")
        os.makedirs(sub, exist_ok=True)
        for d, tag in (("forward", "F"), ("backward", "B")):
            cs = [c for c in self.curves if c.direction == d]
            if not cs:
                continue
            n = max(len(c.V) for c in cs)
            cols, names = [], []
            for c in cs:
                pad = n - len(c.V)
                cols += [np.r_[c.V, np.full(pad, np.nan)], np.r_[c.I, np.full(pad, np.nan)]]
                names += [f"V_{c.index}", f"I_{c.index}"]
            p = os.path.join(sub, f"IV_Curves_{tag}.txt")
            ex.save_columns(p, names, *cols)
            files.append(p)
        m = self.models
        if m is not None:
            p = os.path.join(folder, f"{prefix}_IV_models.txt")
            ex.save_columns(p, ["curve", "forward", "G_low_G0", "Vt_plus_V", "Vt_minus_V", "eps0_TVS_eV", "a_TVS",
                                "eps0_SLM_eV", "eps0_err", "Gamma_SLM_eV", "Gamma_err", "a_SLM", "R2"],
                            m.index, m.direction == "forward", m.g_low, m.v_plus, m.v_minus, m.eps_tvs, m.a_tvs,
                            m.column("eps0"), m.column("eps0_err"), m.column("Gamma"), m.column("Gamma_err"),
                            m.column("a"), m.column("r2"))
            files.append(p)
            if m.most_probable is not None:
                vg, mu, sd, cnt = m.most_probable
                f = m.mean["most probable"][2]
                fit = ivm.slm_current(vg, f.eps0, f.Gamma, f.a, f.n, m.settings.temperature) if f.ok \
                    else np.full(len(vg), np.nan)
                p = os.path.join(folder, f"{prefix}_IV_most_probable.txt")
                ex.save_columns(p, ["V", "I_most_probable_A", "sigma_A", "points", "I_SLM_fit_A"], vg, mu, sd, cnt, fit)
                files.append(p)
            files.append(ex.save_summary(folder, m.summary_lines() + ["settings: " + str(m.settings.to_dict())],
                                         f"{prefix}_IV_models_summary.txt"))
        return files


# Piezo modulation

MOD_LABELS = {"f_mod": "modulation frequency (Hz, 0 = auto)", "periods_per_window": "periods per window",
              "f_search_min": "auto search above (Hz)", "g_min": "log G from", "g_max": "log G to",
              "displacement_ratio": "displacement ratio"}


class ModulationTab(AnalysisTab):
    title = "Piezo modulation"

    def __init__(self, project, main):
        super().__init__(project, main)
        self.rec = QComboBox()
        self.form = SettingsForm("Lock-in demodulation", md.ModulationSettings(), MOD_LABELS)
        self.forms["modulation"] = self.form
        info = QLabel("For recordings with a sinusoidal piezo modulation. ln G and z are demodulated at the "
                      "modulation frequency, and beta = -d(ln G)/dz. Tunnelling gives beta of about 2 kappa "
                      "(about 20 /nm for Au in vacuum, lower in solvents); molecular plateaus give a small beta.")
        info.setWordWrap(True)
        self.result_label = QLabel("")
        self.result_label.setWordWrap(True)
        self.layout_.addWidget(scroll_panel(_recording_box(self, self.rec), self.form,
                                            make_button("Analyse", self.run), self.result_label, info))
        self.plot = PlotPanel()
        self.layout_.addWidget(self.plot, 1)
        self.res = None

    def run(self):
        ms = self.form_value("modulation")
        rec = self.project.recording(self.rec.currentText())
        if ms is None or rec is None:
            return
        self.background(lambda: (rec, md.lockin(rec, ms)), on_done=self._done, busy="Demodulating ...")

    def _done(self, out):
        rec, r = out
        self.res = r
        self.result_label.setText(f"Modulation frequency: {r['f_mod']:.2f} Hz; {len(r['beta'])} windows; "
                                  f"median beta = {np.median(r['beta']) if len(r['beta']) else np.nan:.3g} 1/nm")
        self.rec_used = rec
        self.plot.show_panels(self.panels(), 2)

    def panels(self):
        if self.res is None:
            return []
        return figs.modulation(self.rec_used, self.res, md.spectrum)

    def has_results(self):
        return self.res is not None and len(self.res["beta"]) > 0

    def export(self, folder, prefix):
        if not self.has_results():
            return []
        r = self.res
        f = os.path.join(folder, f"{prefix}_piezo_modulation.txt")
        ex.save_columns(f, ["t_s", "logG", "amp_z_nm", "amp_lnG", "beta_per_nm", "phase_deg"],
                        r["t"], r["logG"], r["amp_z"], r["amp_lnG"], r["beta"], r["phase"])
        return [f]
