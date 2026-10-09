"""Start script for the GUI.

Run it directly, or use the exe made by build_exe.py.
"BJ_Analysis.exe --selftest report.txt" runs everything on simulated data
without opening a window and writes a short report. I use it to check a
fresh build.
"""

import os
import sys


def selftest(report_path):
    import numpy as np
    from break_junction import analysis as an, clustering as cl, noise as nz, simulation as sim, traces as tr
    from break_junction import G0, iv
    from break_junction.conversion import Recording

    lines = []
    rec = sim.simulate(sim.SimulationSettings(n_cycles=120, second_molecule_logG=-2.2, p_second=0.4))
    T, H = tr.extract_traces(rec, tr.SegmentSettings())
    good = tr.selected(T)
    lines.append(f"traces {len(T)} selected {len(good)} holds {len(H)}")
    x, y = an.hist1d(good, an.HistSettings())
    peaks, *_ = an.fit_peaks(x, y, -4.5, -1.5, 2)
    lines.append("peaks " + ", ".join(f"{p.center:.2f}" for p in peaks))
    res = cl.run(good, cl.ClusterSettings(n_clusters=3), scan_k=True)
    lines.append(f"clusters {np.bincount(res.labels).tolist()} silhouette {res.silhouette:.3f}")
    ns = nz.morris2025_settings(nz.NoiseSettings(g_min=-5.0, g_max=-1.0))
    r = nz.analyse(H, ns)
    G, NP = r.arrays()
    sc = nz.scaling_exponent(G, NP, ns)
    lines.append(f"noise accepted {len(G)} n_TSE {sc['n_tse']:.3f} n_OLS {sc['n_ols']:.3f} (truth 1.5)")
    t = np.arange(20000) / 1e4
    V = 0.5 * (2 / np.pi) * np.arcsin(np.sin(2 * np.pi * 2 * t))
    I = G0 * 1e-3 * V
    curves = iv.curves_from_recording(Recording("iv", 1e4, np.full(len(t), 1e-3), None, V, I), iv.IVSettings())
    lines.append(f"iv curves {len(curves)}")
    from break_junction import bj_io, conversion as cv
    logG = np.linspace(0.5, -6, 5000)
    Vb = 0.1 + 0.01 * np.sin(np.arange(5000) / 50)
    I = 10 ** logG * G0 * Vb
    table = bj_io.ColumnTable("dual", {"Ia": np.clip(I / 1e-5, -10.5, 10.5), "Ib": np.clip(I / 1e-9, -10, 10),
                                       "V": Vb})
    rec = cv.table_to_recordings(table, cv.ColumnSettings(
        amplifier="dual stage", current_a_column="Ia", current_b_column="Ib", current_scale_a=1e-5,
        current_scale_b=1e-9, b_saturation=9.99, voltage_column="V", piezo_column=""))[0]
    err = np.max(np.abs(np.log10(rec.G) - logG))
    import igor2.packed  # noqa: F401 - must be bundled for .pxp files
    lines.append(f"dual-stage columns max error {err:.1e} decades, igor2 ok")
    import tempfile
    from break_junction import figures, style
    zc, gc, H2 = an.hist2d(good, an.HistSettings())
    panels = [figures.histogram_1d(x, y, "per trace", "All accepted", len(good)),
              figures.histogram_2d(zc, gc, H2, "per trace")] + figures.noise(r, sc, ns, H[0], 0)
    with tempfile.TemporaryDirectory() as tmp:
        files = style.export_panels(panels, tmp, "selftest", style.FigureSettings())
        sizes = [os.path.getsize(f) for f in files]
    lines.append(f"figures {len(files)} files, all non-empty: {all(s > 1000 for s in sizes)}")
    lines.append("SELFTEST OK")
    with open(report_path, "w") as fh:
        fh.write("\n".join(lines) + "\n")
    return 0


if __name__ == "__main__":
    # make the package importable no matter where this file is started from
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    if len(sys.argv) >= 3 and sys.argv[1] == "--selftest":
        try:
            code = selftest(sys.argv[2])
        except Exception:  # noqa: BLE001 - the built program has no console: report in the file
            import traceback
            with open(sys.argv[2], "w") as fh:
                fh.write("SELFTEST FAILED\n" + traceback.format_exc())
            code = 1
        sys.exit(code)
    from break_junction.gui.app import main

    sys.exit(main())
