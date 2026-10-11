"""Start script for the GUI.

Run it directly, or use the exe made by build_exe.py.
"GZero.exe --selftest report.txt" runs everything on simulated data
without opening a window and writes a short report. I use it to check a
fresh build.
"""

import os
import sys


def selftest(report_path):
    import numpy as np
    from gzero import analysis as an, clustering as cl, noise as nz, simulation as sim, traces as tr
    from gzero import G0, iv
    from gzero.conversion import Recording

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
    from gzero import ivmodels as ivm, junctions as jn
    jr = jn.analyse(good, jn.JunctionSettings())
    lines.append(f"junction yield {100 * jr.yield_()[0]:.1f} % (70 % molecular; the -2.2 class sits at the window edge)")
    t = np.arange(20000) / 1e4
    V = 1.5 * (2 / np.pi) * np.arcsin(np.sin(2 * np.pi * 2 * t))
    I = ivm.slm_current(V, 0.8, 0.02)
    curves = iv.curves_from_recording(Recording("iv", 1e4, np.full(len(t), 1e-3), None, V, I), iv.IVSettings())
    lines.append(f"iv curves {len(curves)}")
    mr = ivm.analyse(curves, ivm.IVModelSettings(), iv.mean_curve)
    f = mr.mean["forward"][2]
    lines.append(f"SLM eps0 {f.eps0:.3f} eV Gamma {1e3 * f.Gamma:.1f} meV (truth 0.800, 20.0); "
                 f"eps0 from V_t {np.nanmedian(mr.eps_tvs):.3f} eV")
    s3 = ivm.rashid2025_settings()
    Vg = np.linspace(-1, 1, 81)
    f3 = ivm.fit_slm(Vg, ivm.slm_current(Vg, 1.10, 0.0113, 0, 1, 300), s3, grid=True)
    lines.append(f"SLM 300 K energy grid + LM: eps0 {f3.eps0:.3f} eV Gamma {1e3 * f3.Gamma:.2f} meV (truth 1.100, 11.30)")
    dw = an.dwell_time_ms(peaks[0], 1 / an.HistSettings().bins_per_decade, rec.fs)
    lines.append(f"dwell time of the first peak {dw:.1f} ms per trace")
    from gzero import egain as eg, events as evm
    er = evm.analyse(H, evm.EventSettings(source="holds"))
    lines.append(f"events in holds: {int(er.n_flicker.sum())} flickers, {int(er.n_mechanical.sum())} mechanical")
    sets = [eg.analyse(eg.Dataset(f"C{L}", eg.simulate(length=L, seed=L), L), eg.EGaInSettings()) for L in (10, 14, 18)]
    b = eg.fit_beta(sets, -0.5)
    lines.append(f"EGaIn beta {b['beta']:.3f} +/- {b['beta_err']:.3f} (truth 0.9), yield {sets[0].n_working}/{sets[0].n_junctions}")
    from gzero import thermo as th
    tr_ = th.analyse(th.ThermoDataset("sim", th.simulate(S=10.0, seed=0)), th.ThermoSettings(bootstrap=20))
    lines.append(f"thermopower S {tr_.S:.2f} +/- {tr_.S_err:.2f} uV/K (truth 10)")
    import tempfile as _tf
    with _tf.TemporaryDirectory() as tmp:
        js = eg.simulate(n_junctions=3, length=12, seed=1, p_short=0)
        cols = [c for j in js for c in (j.V, j.J * eg.EGaInSettings().area())]
        np.savetxt(os.path.join(tmp, "pairs.txt"), np.column_stack(cols), delimiter="\t")
        jp = eg.load_junctions(os.path.join(tmp, "pairs.txt"), eg.EGaInSettings(layout=eg.LAYOUTS[1]))
    lines.append(f"EGaIn column pairs: {len(jp)} junctions from one file")
    # a big file read in pieces must give the same traces as the whole file
    from nptdms import TdmsWriter, ChannelObject
    from gzero import streaming as stm, conversion as cv0
    srec = sim.simulate(sim.SimulationSettings(n_cycles=30, seed=7))
    with _tf.TemporaryDirectory() as tmp:
        tp = os.path.join(tmp, "s.tdms")
        with TdmsWriter(tp) as wr:
            props = {"wf_increment": 1.0 / srec.fs}
            for a0 in range(0, srec.n, 200_000):
                sl = slice(a0, a0 + 200_000)
                wr.write_segment([ChannelObject("BJRawData", "M_E1", np.full(len(srec.G[sl]), 0.1), properties=props),
                                  ChannelObject("BJRawData", "M_i1a", srec.G[sl] * G0 * 0.1 / 1e-5,
                                                properties={**props, "scale_factor": 1e-5}),
                                  ChannelObject("BJRawData", "U_Piezo", srec.piezo[sl], properties=props)])
        ccs = cv0.ConversionSettings(amplifier="single stage")
        full = cv0.load_recordings(tp, None, ccs)[0]
        tf_, hf_ = tr.extract_traces(full, tr.SegmentSettings())
        light = stm.load_streamed(tp, ccs, stm.StreamSettings(chunk_s=15, margin_s=8, overview_points=5000))
        ts_, hs_ = stm.stream_traces(light, tr.SegmentSettings())
        same = [(x.start, x.stop) for x in tf_] == [(x.start, x.stop) for x in ts_] and len(hf_) == len(hs_)
    lines.append(f"TDMS read in pieces: {len(ts_)} traces, {len(hs_)} holds, same as the whole file: {same}")
    if not same:
        raise AssertionError("reading in pieces changed the traces")
    from gzero import bj_io, conversion as cv
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
    from gzero import figures, style
    zc, gc, H2 = an.hist2d(good, an.HistSettings())
    panels = [figures.histogram_1d(x, y, "per trace", "All accepted", len(good)),
              figures.histogram_2d(zc, gc, H2, "per trace")] + figures.noise(r, sc, ns, H[0], 0) \
        + figures.junctions(jr) + figures.iv_models(mr, ivm.slm_current, ivm.fowler_nordheim) \
        + figures.events(er, None, evm.time_histogram(H, evm.EventSettings())) \
        + figures.egain(sets[0], b, ivm.slm_current) + figures.thermo(tr_)
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
    from gzero.gui.app import main

    sys.exit(main())
