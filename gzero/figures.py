"""Plots used both on screen and for the figure export.

Each function returns Panel objects. The GUI puts them in a grid, the
export saves every one as a separate figure (style.export_panels).
"""

from __future__ import annotations

import numpy as np

from .style import PALETTE, Panel, is_export, title

LOGG = "log(G/G$_0$)"


def decimate(x: np.ndarray, y: np.ndarray, max_points: int = 200_000):
    """Min/max decimation so long recordings plot fast without losing spikes."""
    n = len(y)
    if n <= max_points:
        return x, y
    k = int(np.ceil(n / (max_points / 2)))
    m = n // k
    yr = y[:m * k].reshape(m, k)
    xr = x[:m * k].reshape(m, k)
    imin, imax = yr.argmin(axis=1), yr.argmax(axis=1)
    rows = np.arange(m)
    first = np.minimum(imin, imax)
    second = np.maximum(imin, imax)
    xs = np.column_stack([xr[rows, first], xr[rows, second]]).ravel()
    ys = np.column_stack([yr[rows, first], yr[rows, second]]).ravel()
    return xs, ys


def _color(k: int) -> str:
    return PALETTE[k % len(PALETTE)]


def _norm_label(normalize: str) -> str:
    return {"per trace": "Counts per trace", "density": "Probability density"}.get(normalize, "Counts")


def _colorbar(fig, ax, mappable, label):
    cb = fig.colorbar(mappable, ax=ax, label=label, fraction=0.05, pad=0.02)
    cb.outline.set_linewidth(ax.spines["left"].get_linewidth())
    return cb


# Traces and histograms

def trace_overlay(traces, n_max: int = 50, seed: int = 0) -> list[Panel]:
    def draw(fig, ax):
        if not traces:
            return
        rng = np.random.default_rng(seed)
        pick = rng.choice(len(traces), min(n_max, len(traces)), replace=False)
        for i in pick:
            ax.plot(traces[i].z, traces[i].logG, lw=0.4 if is_export(ax) else 0.5, alpha=0.6)
        ax.set_xlabel("Distance (nm)")
        ax.set_ylabel(LOGG)
        title(ax, f"{len(pick)} random accepted traces")
    return [Panel("example_traces", draw)]


def histogram_1d(x, y, normalize, label, n, peaks=None, clusters=None) -> Panel:
    """1D log-G histogram; peaks = (PeakFit list, xf, yf, base); clusters = [(label, x, y)]."""
    def draw(fig, ax):
        ax.fill_between(x, y, step="mid", color="0.75", lw=0)
        ax.step(x, y, where="mid", color="k", lw=0.6, label=f"{label} ({n} traces)")
        for k, (cl, xk, yk) in enumerate(clusters or []):
            ax.step(xk, yk, where="mid", lw=0.8, color=_color(k), label=cl)
        if peaks:
            pk, xf, yf, _ = peaks
            ax.plot(xf, yf, color=PALETTE[3], lw=1.0, label="Gaussian fit")
            for p in pk:
                ax.axvline(p.center, color=PALETTE[3], ls="--", lw=0.5)
                ax.annotate(f"{p.G:.2e} G$_0$", (p.center, np.interp(p.center, xf, yf)),
                            textcoords="offset points", xytext=(3, 4), color=PALETTE[3],
                            fontsize=None if is_export(ax) else 8)
        ax.set_xlabel(LOGG)
        ax.set_ylabel(_norm_label(normalize))
        ax.set_ylim(bottom=0)
        ax.legend(loc="best")
    return Panel("logG_histogram", draw)


def histogram_2d(zc, gc, H, normalize, cmap="viridis", vmax_pct=99.5, name="2D_histogram") -> Panel:
    def draw(fig, ax):
        vmax = np.percentile(H[H > 0], vmax_pct) if np.any(H > 0) else 1
        m = ax.pcolormesh(zc, gc, H.T, cmap=cmap, vmin=0, vmax=vmax, shading="nearest", rasterized=True)
        _colorbar(fig, ax, m, _norm_label(normalize))
        ax.set_xlabel("Distance (nm)")
        ax.set_ylabel(LOGG)
        title(ax, "2D conductance-distance histogram")
    return Panel(name, draw, aspect=0.8)


# Plateau length

def plateau(windows, lengths, max_length, bins, min_length, betas=None, beta_median=None) -> list[Panel]:
    edges = np.linspace(0, max_length, bins + 1)

    def draw_len(fig, ax):
        for k, ((lo, hi), L) in enumerate(zip(windows, lengths)):
            Lk = L[L > min_length]
            mean = np.mean(Lk) if len(Lk) else np.nan
            ax.hist(Lk, edges, histtype="step", lw=1.0, color=_color(k),
                    label=f"[{lo:g}, {hi:g}]: mean {mean:.3f} nm")
        ax.set_xlabel("Plateau length (nm)")
        ax.set_ylabel("Traces")
        ax.legend(title=f"{LOGG} window")

    panels = [Panel("plateau_length_histogram", draw_len)]
    if betas is not None and len(betas):
        def draw_beta(fig, ax):
            ax.hist(betas, 40, color="0.6")
            ax.axvline(beta_median, color=PALETTE[3], lw=1.0, label=f"median {beta_median:.3g} nm$^{{-1}}$")
            ax.set_xlabel(r"Tunnelling decay $\beta$ (nm$^{-1}$)")
            ax.set_ylabel("Traces")
            ax.legend()
        panels.append(Panel("tunnelling_decay_histogram", draw_beta))
    return panels


# Correlation

def correlation_map(c, C, limit, cmap) -> Panel:
    def draw(fig, ax):
        m = ax.pcolormesh(c, c, C, cmap=cmap, vmin=-limit, vmax=limit, shading="nearest", rasterized=True)
        _colorbar(fig, ax, m, "Correlation")
        ax.set_xlabel(LOGG)
        ax.set_ylabel(LOGG)
        ax.set_aspect("equal")
    return Panel("correlation_map", draw, aspect=0.85)


def mean_histogram(c, h) -> Panel:
    def draw(fig, ax):
        ax.step(c, h, where="mid", color="k")
        ax.set_xlabel(LOGG)
        ax.set_ylabel("Counts per trace")
        ax.set_ylim(bottom=0)
        title(ax, "Mean per-trace histogram")
    return Panel("correlation_mean_histogram", draw, aspect=0.6)


# Clustering

def clustering(res, cluster_hists, cluster_2d) -> list[Panel]:
    """cluster_hists = [(x, y, n)], cluster_2d = [(zc, gc, H, n)] per cluster."""
    k = res.n_clusters

    def scree(fig, ax):
        n = len(res.explained)
        ax.bar(np.arange(1, n + 1), 100 * res.explained, color="0.6")
        ax.plot(np.arange(1, n + 1), 100 * np.cumsum(res.explained), "k.-", lw=0.6)
        ax.set_xlabel("Principal component")
        ax.set_ylabel("Explained variance (%)")
        title(ax, "Scree plot")

    def scatter(fig, ax):
        for j in range(k):
            sel = res.labels == j
            ax.scatter(res.scores[sel, 0], res.scores[sel, 1], s=2 if is_export(ax) else 4,
                       color=_color(j), label=f"Cluster {j + 1}", linewidths=0)
        ax.set_xlabel("PC 1")
        ax.set_ylabel("PC 2")
        ax.legend(markerscale=3)

    def hists(fig, ax):
        for j, (x, y, n) in enumerate(cluster_hists):
            ax.plot(x, y, color=_color(j), lw=0.9, label=f"Cluster {j + 1} ({n})")
        ax.set_xlabel(LOGG)
        ax.set_ylabel("Counts per trace")
        ax.set_ylim(bottom=0)
        ax.legend()

    panels = [Panel("clustering_PCA_scree", scree), Panel("clustering_PCA_scatter", scatter), Panel("clustering_logG_histograms", hists)]
    if res.k_scan:
        def scan(fig, ax):
            ks = sorted(res.k_scan)
            ax.plot(ks, [res.k_scan[q] for q in ks], "o-", color="k", ms=3)
            ax.set_xlabel("Number of clusters k")
            ax.set_ylabel("Silhouette score")
        panels.append(Panel("clustering_silhouette_scan", scan))
    for j, (zc, gc, H, n) in enumerate(cluster_2d):
        def draw2d(fig, ax, zc=zc, gc=gc, H=H, n=n, j=j):
            vmax = np.percentile(H[H > 0], 99.5) if np.any(H > 0) else 1
            m = ax.pcolormesh(zc, gc, H.T, cmap="viridis", vmin=0, vmax=vmax, shading="nearest", rasterized=True)
            _colorbar(fig, ax, m, "Counts per trace")
            ax.set_xlabel("Distance (nm)")
            ax.set_ylabel(LOGG)
            title(ax, f"Cluster {j + 1}: {n} traces", color=_color(j))
        panels.append(Panel(f"clustering_cluster{j + 1}_2D_histogram", draw2d, aspect=0.8))
    return panels


# Flicker noise

def noise(res, sc, ns, segment=None, segment_index=0) -> list[Panel]:
    panels = []
    if segment is not None:
        def seg(fig, ax):
            t, y = decimate(segment.t, segment.logG)
            ax.plot(t, y, lw=0.4, color="k")
            for w in res.windows:
                if w.segment == segment_index:
                    ax.axvspan(w.start / segment.fs, (w.start + w.length) / segment.fs,
                               color=PALETTE[2] if w.accepted else PALETTE[3], alpha=0.15, lw=0)
            ax.set_xlabel("Time (s)")
            ax.set_ylabel(LOGG)
            title(ax, f"Segment {segment_index}: green = accepted, red = rejected")
        panels.append(Panel("noise_example_segment", seg, aspect=0.6))

    if len(res.freqs):
        def psd(fig, ax):
            f = res.freqs[1:]
            ax.loglog(f, res.psd_mean[1:], color="k", lw=0.6, label="Mean of accepted windows")
            if res.psd_floor is not None:
                ax.loglog(f, res.psd_floor[1:], color="0.6", lw=0.6, label="Noise floor")
            ax.axvspan(ns.f_lo, ns.f_hi, color=PALETTE[5], alpha=0.15, lw=0)
            band = (f >= ns.f_lo) & (f <= ns.f_hi)
            if np.isfinite(res.alpha) and np.any(band):
                f0, p0 = f[band][0], res.psd_mean[1:][band][0]
                ax.loglog(f[band], p0 * (f[band] / f0) ** (-res.alpha), "--", color=PALETTE[3], lw=0.8,
                          label=f"1/f$^{{{res.alpha:.2f}}}$")
            ax.set_xlabel("Frequency (Hz)")
            ax.set_ylabel("PSD (G$_0^2$ Hz$^{-1}$)")
            ax.legend()
        panels.append(Panel("noise_PSD", psd))

    if not sc:
        return panels
    x, ly = sc["x"], sc["log_np"]

    def fits(fig, ax):
        ax.plot(x, ly, ".", ms=1.5 if is_export(ax) else 2, color=PALETTE[0], alpha=0.5)
        xx = np.array([x.min(), x.max()])
        ax.plot(xx, sc["ols_intercept"] + sc["n_ols"] * xx, color=PALETTE[1], lw=1.0,
                label=f"OLS, n = {sc['n_ols']:.2f} $\\pm$ {sc['n_ols_se']:.2f}")
        ax.plot(xx, sc["tse_intercept"] + sc["n_tse"] * xx, color="k", lw=1.0,
                label=f"Theil-Sen, n = {sc['n_tse']:.2f} $\\pm$ {sc['n_tse_se']:.2f}")
        ax.set_xlabel(LOGG)
        ax.set_ylabel("log(NP/G$_0^2$)")
        ax.legend()

    def bell(fig, ax):
        from .noise import _gauss2d
        H, xe, ye, b = sc["hist"], sc["x_edges"], sc["y_edges"], sc["bell"]
        m = ax.pcolormesh(0.5 * (xe[1:] + xe[:-1]), 0.5 * (ye[1:] + ye[:-1]), H.T, cmap="viridis",
                          shading="nearest", rasterized=True)
        _colorbar(fig, ax, m, "Counts")
        X, Y = np.meshgrid(np.linspace(xe[0], xe[-1], 120), np.linspace(ye[0], ye[-1], 120))
        Z = _gauss2d((X, Y), b.amplitude, b.x0, b.y0, b.sx, b.sy, b.rho, 0).reshape(X.shape)
        ax.contour(X, Y, Z, levels=b.amplitude * np.array([0.2, 0.5, 0.8]), colors="w",
                   linewidths=0.5 if is_export(ax) else 0.8)
        ax.set_xlabel(LOGG)
        ax.set_ylabel(f"log(NP/G$^{{{sc['n_best']:.2f}}}$)")
        title(ax, f"Bell at n = {sc['n_best']:.2f} (fitted correlation {b.rho:.2f})")

    def corr(fig, ax):
        ax.plot(sc["n_grid"], sc["rho_pearson"], color=PALETTE[1], label="Pearson r (OLS)")
        ax.plot(sc["n_grid"], sc["tau_kendall"], color="k", label=r"Kendall $\tau$ (Theil-Sen)")
        ax.plot(sc["n_grid"], sc["rho_fit"], color=PALETTE[3], ls="--", label="2D Gaussian fit")
        ax.axhline(0, color="0.6", lw=0.5)
        ax.axvline(sc["n_tse"], color="k", ls=":", lw=0.6)
        ax.set_xlabel("Scaling exponent n")
        ax.set_ylabel("Correlation of log G and log(NP/G$^n$)")
        ax.legend()

    def resid(fig, ax):
        r_ols = ly - (sc["ols_intercept"] + sc["n_ols"] * x)
        r_tse = ly - (sc["tse_intercept"] + sc["n_tse"] * x)
        parts = ax.violinplot([r_ols, r_tse], showmedians=True)
        for body, c in zip(parts["bodies"], (PALETTE[1], "0.4")):
            body.set_facecolor(c)
        for key in ("cbars", "cmins", "cmaxes", "cmedians"):
            if key in parts:
                parts[key].set_color("k")
                parts[key].set_linewidth(0.6)
        ax.set_xticks([1, 2], ["OLS", "Theil-Sen"])
        ax.axhline(0, color="0.6", lw=0.5)
        ax.set_ylabel("Residual of log(NP)")
        title(ax, "Fit residuals (long tails indicate outliers)")

    panels += [Panel("noise_scaling_fits", fits), Panel("noise_bell_histogram", bell, aspect=0.8),
               Panel("noise_correlation_vs_n", corr), Panel("noise_fit_residuals", resid)]
    return panels


# I-V and piezo modulation

def iv(curves, s, hist_fn, mean_fn) -> list[Panel]:
    def overlay(fig, ax):
        for c in curves[:200]:
            ax.plot(c.V, c.I, lw=0.3, alpha=0.5, color=PALETTE[0] if c.direction == "forward" else PALETTE[1])
        for d, (vg, m, _) in mean_fn(curves).items():
            ax.plot(vg, m, color="k" if d == "forward" else PALETTE[3], lw=1.0, label=f"Mean {d}")
        ax.set_xlabel("Bias (V)")
        ax.set_ylabel("Current (A)")
        ax.legend()

    panels = [Panel("IV_curves", overlay)]
    for q, name, lab in (("I", "IV_current_histogram", "log|I| (A)" if s.log_current else "Current (A)"),
                         ("G", "IV_conductance_histogram", "log(G/G$_0$), G = I/V"),
                         ("dI/dV", "IV_dIdV_histogram", "log[(dI/dV)/G$_0$]")):
        def draw(fig, ax, q=q, lab=lab):
            H, ve, ye = hist_fn(curves, s, q)
            m = ax.pcolormesh(0.5 * (ve[1:] + ve[:-1]), 0.5 * (ye[1:] + ye[:-1]), H.T, cmap="viridis",
                              shading="nearest", vmax=np.percentile(H[H > 0], 99) if np.any(H > 0) else 1,
                              rasterized=True)
            _colorbar(fig, ax, m, "Counts per curve")
            ax.set_xlabel("Bias (V)")
            ax.set_ylabel(lab)
        panels.append(Panel(name, draw, aspect=0.8))
    return panels


def modulation(rec, r, spectrum_fn) -> list[Panel]:
    n = min(rec.n, int(rec.fs * 4))

    def spec_piezo(fig, ax):
        f, a = spectrum_fn(rec.piezo[:n], rec.fs)
        ax.semilogy(f[1:], a[1:], lw=0.5, color=PALETTE[4])
        ax.axvline(r["f_mod"], color=PALETTE[3], lw=0.5)
        ax.set_xlabel("Frequency (Hz)")
        ax.set_ylabel("Piezo amplitude (nm)")

    def spec_g(fig, ax):
        f, a = spectrum_fn(rec.logG()[:n], rec.fs)
        ax.semilogy(f[1:], a[1:], lw=0.5, color="k")
        ax.axvline(r["f_mod"], color=PALETTE[3], lw=0.5)
        ax.set_xlabel("Frequency (Hz)")
        ax.set_ylabel("Amplitude of log G")

    panels = [Panel("modulation_piezo_spectrum", spec_piezo), Panel("modulation_logG_spectrum", spec_g)]
    if len(r["beta"]):
        def beta_g(fig, ax):
            ax.plot(r["logG"], r["beta"], ".", ms=1.5, alpha=0.5, color=PALETTE[0])
            ax.set_xlabel(LOGG)
            ax.set_ylabel(r"$\beta$ (nm$^{-1}$)")

        def beta_h(fig, ax):
            ax.hist(r["beta"], 60, color="0.6")
            ax.set_xlabel(r"$\beta$ (nm$^{-1}$)")
            ax.set_ylabel("Windows")
        panels += [Panel("modulation_beta_vs_G", beta_g), Panel("modulation_beta_histogram", beta_h)]
    return panels


def iv_models(res, slm_fn, fn_fn) -> list[Panel]:
    """Single-level fit and transition voltage spectroscopy (ivmodels.analyse result)."""
    panels = []
    colors = {"forward": "k", "backward": PALETTE[3]}
    if res.mean:
        def slm(fig, ax):
            for d, (vg, m, f, _, _) in res.mean.items():
                ax.plot(vg, m * 1e9, ".", ms=1.5 if is_export(ax) else 2.5, color=colors[d], alpha=0.5,
                        label=f"Mean {d}")
                if f.ok:
                    ax.plot(vg, slm_fn(vg, f.eps0, f.Gamma, f.a, res.settings.n_molecules) * 1e9,
                            color=PALETTE[0] if d == "forward" else PALETTE[1], lw=1.0,
                            label=f"SLM: $\\varepsilon_0$ = {f.eps0:.2f} eV, $\\Gamma$ = {1e3 * f.Gamma:.2g} meV")
            ax.set_xlabel("Bias (V)")
            ax.set_ylabel("Current (nA)")
            ax.legend()

        def fn(fig, ax):
            for d, (vg, m, _, vp, vm) in res.mean.items():
                for side, (x, y) in fn_fn(vg, m, res.settings.tvs_v_min).items():
                    ax.plot(x, y, lw=0.8, color=colors[d], ls="-" if side == "positive" else "--",
                            label=f"{d}, {side} bias")
                for vt in (vp, vm):
                    if np.isfinite(vt):
                        ax.axvline(1 / vt, color=colors[d], lw=0.5, ls=":")
            ax.set_xlabel("1/V (V$^{-1}$)")
            ax.set_ylabel("ln(|I|/V$^2$)")
            ax.legend(loc="lower right")
            title(ax, "Fowler-Nordheim plot (minimum at 1/V$_t$)")
        panels += [Panel("IV_SLM_fit", slm), Panel("IV_Fowler_Nordheim", fn)]

    vp, vm = res.v_plus[np.isfinite(res.v_plus)], res.v_minus[np.isfinite(res.v_minus)]
    if len(vp) or len(vm):
        def vt(fig, ax):
            allv = np.r_[vp, np.abs(vm)]
            edges = np.linspace(allv.min(), allv.max(), 40) if np.ptp(allv) > 0 else 20
            if len(vp):
                ax.hist(vp, edges, histtype="step", color=PALETTE[0], lw=1.0, label=f"V$_t^+$, median {np.median(vp):.2f} V")
            if len(vm):
                ax.hist(np.abs(vm), edges, histtype="step", color=PALETTE[1], lw=1.0,
                        label=f"|V$_t^-$|, median {np.median(np.abs(vm)):.2f} V")
            ax.set_xlabel("Transition voltage (V)")
            ax.set_ylabel("Curves")
            ax.legend()
        panels.append(Panel("IV_transition_voltage_histogram", vt))

    eps, gam = res.column("eps0"), res.column("Gamma")
    ok = np.isfinite(eps) & np.isfinite(gam) & (gam > 0)
    if ok.sum() >= 2:
        def params(fig, ax):
            ax.scatter(eps[ok], 1e3 * gam[ok], s=2 if is_export(ax) else 5, color=PALETTE[0], linewidths=0,
                       alpha=0.7)
            ax.set_yscale("log")
            ax.set_xlabel(r"$\varepsilon_0$ (eV)")
            ax.set_ylabel(r"$\Gamma$ (meV)")
            title(ax, f"SLM fit of {ok.sum()} curves")

        def eps_hist(fig, ax):
            et = res.eps_tvs[np.isfinite(res.eps_tvs)]
            both = np.r_[eps[ok], et]
            edges = np.linspace(both.min(), both.max(), 40) if np.ptp(both) > 0 else 20
            ax.hist(eps[ok], edges, histtype="step", color=PALETTE[0], lw=1.0,
                    label=f"SLM fit, median {np.median(eps[ok]):.2f} eV")
            if len(et):
                ax.hist(et, edges, histtype="step", color=PALETTE[2], lw=1.0,
                        label=f"from V$_t$, median {np.median(et):.2f} eV")
            ax.set_xlabel(r"$\varepsilon_0$ (eV)")
            ax.set_ylabel("Curves")
            ax.legend()
        panels += [Panel("IV_SLM_parameters", params), Panel("IV_eps0_histogram", eps_hist)]
    return panels


# Junction statistics

def junctions(res) -> list[Panel]:
    s = res.settings
    j = res.junction
    p, se = res.yield_()
    bl = res.blocks()

    def g_hist(fig, ax):
        g = res.logG[j]
        ax.hist(g, np.linspace(s.g_lo, s.g_hi, s.g_bins + 1), color="0.6")
        if len(g):
            ax.axvline(np.median(g), color=PALETTE[3], lw=1.0,
                       label=f"median {np.median(g):.2f} ({10 ** np.median(g):.2e} G$_0$)")
            ax.legend()
        ax.set_xlabel(f"Plateau {LOGG}")
        ax.set_ylabel("Traces")
        title(ax, f"{int(j.sum())} of {res.n} traces with a plateau")

    def g_vs_n(fig, ax):
        ax.plot(res.index[j], res.logG[j], ".", ms=1.5 if is_export(ax) else 2.5, color=PALETTE[0], alpha=0.5)
        if len(bl):
            ax.plot(bl[:, 0], bl[:, 3], "o-", color="k", ms=3, lw=0.8, label=f"Median per {s.block_size} traces")
            ax.legend()
        ax.set_xlabel("Trace number")
        ax.set_ylabel(f"Plateau {LOGG}")

    def yield_vs_n(fig, ax):
        if len(bl):
            ax.errorbar(bl[:, 0], 100 * bl[:, 1], yerr=100 * bl[:, 2], fmt="o-", color="k", ms=3, lw=0.8,
                        capsize=2)
        ax.axhline(100 * p, color=PALETTE[3], lw=0.8, ls="--", label=f"All: {100 * p:.1f} $\\pm$ {100 * se:.1f} %")
        ax.set_ylim(0, 100)
        ax.set_xlabel("Trace number")
        ax.set_ylabel("Junction yield (%)")
        ax.legend()

    def length_hist(fig, ax):
        top = max(1.0, float(np.percentile(res.length, 99))) if res.n else 1.0
        ax.hist(res.length, np.linspace(0, top, 50), color="0.6")
        ax.axvline(s.min_plateau, color=PALETTE[3], lw=0.8, ls="--", label=f"threshold {s.min_plateau:g} nm")
        ax.set_xlabel("Flat plateau length (nm)")
        ax.set_ylabel("Traces")
        ax.legend()

    def length_vs_g(fig, ax):
        ax.plot(res.logG[j], res.length[j], ".", ms=1.5 if is_export(ax) else 2.5, color=PALETTE[0], alpha=0.6)
        ax.set_xlabel(f"Plateau {LOGG}")
        ax.set_ylabel("Plateau length (nm)")

    return [Panel("junction_plateau_conductance_histogram", g_hist), Panel("junction_conductance_vs_trace", g_vs_n),
            Panel("junction_yield_vs_trace", yield_vs_n), Panel("junction_plateau_length_histogram", length_hist),
            Panel("junction_length_vs_conductance", length_vs_g)]
