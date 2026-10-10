# GZero

**Break-junction analysis software.** GZero (G-zero, after the conductance quantum
G0 = 2e²/h) analyses single-molecule break-junction data (MCBJ and STM-BJ): raw TDMS, text
(.cvr/.txt/.csv) or Igor (.pxp/.ibw) files to conductance, trace detection, histograms, plateau
length, junction yield, dwell time, 2D correlation, clustering, flicker noise, flickering and mechanical
events, I-V (single-level model fits at finite temperature, transition voltage spectroscopy), piezo
modulation, and ensemble statistics of large-area (EGaIn) junctions.

The manual is `docs/User_Manual.pdf` (Help -> User manual in the program).

## Getting it

**Windows, no Python needed:** download `GZero-v1.2-windows-x64.zip` from the
[latest release](https://github.com/arman-rashid/GZero/releases/latest), unzip it and
double-click `GZero\GZero.exe`. Keep the folder together.

**Source code:**

```
git clone https://github.com/arman-rashid/GZero.git
cd GZero
```

or download the ZIP from the green Code button on GitHub.

## Running it

From Python (3.9+), in the project folder (the one with `GZero.py`):

```
py -m pip install -r requirements.txt
py GZero.py
```

(`py` is the Windows launcher; on macOS / Linux use `python3`.)

Windows exe, no Python needed afterwards:

```
py -m pip install pyinstaller
py build_exe.py
```

The program ends up in `dist/GZero/GZero.exe`. Copy the whole `GZero` folder,
the exe doesn't work on its own.

`GZero.exe --selftest report.txt` runs everything on simulated data and writes a short
report. Last line should be `SELFTEST OK`.

File -> Save results writes the results as text files and every plot as a separate figure
(Nature format: 89 mm wide, Arial 7 pt, pdf + 600 dpi png) into a `figures` subfolder. The
figure settings are under File -> Figure export settings. There "plot with" can be set to Origin:
then every plot is rebuilt in Origin instead (data in workbooks, graphs with the same size and
style, saved as `<name>.opju`). File -> Send plots to Origin does this without saving files. This
needs Windows, Origin 2021+ and `py -m pip install originpro`.

## Project layout

```
GZero.py      start the program (also: py -m gzero)
build_exe.py        build the Windows exe with PyInstaller
requirements.txt
gzero/     the Python package: analysis modules and gui/
docs/               user manual (Markdown source, PDF, build_manual.py)
```

## Code

The modules in `gzero/` work without the GUI too, e.g.

```python
from gzero import conversion as cv, traces as tr
rec = cv.load_recording("measurement.tdms")
traces, holds = tr.extract_traces(rec, tr.SegmentSettings())
```

Modules:

| Module | |
|---|---|
| `bj_io` | reading TDMS / .cvr / text / Igor .pxp and .ibw, writing .cvr |
| `conversion` | calibration (linear or diode polynomial), single- or dual-stage amplifier (merging Ia and Ib), series resistor correction, filters; column / channel mapping, so conductance, voltage, piezo, current and time given in a file are used as they are |
| `traces` | splitting recordings into traces, holds, cutting limits, distance axis |
| `analysis` | histograms, Gaussian fits, plateau length, correlation map, tunnelling decay |
| `junctions` | per-trace plateau detection (slope test), junction yield, plateau conductance per trace, drift over the run |
| `origin` | rebuilding the plots in Origin (OriginLab) from the drawn matplotlib panels |
| `clustering` | 2D-histogram features, PCA, k-means / Gaussian mixture / Ward / spectral |
| `noise` | noise power, window selection, scaling exponent n (Theil-Sen, OLS, 2D Gaussian) |
| `stationarity` | ADF test, same results as statsmodels' adfuller |
| `iv`, `modulation`, `simulation` | I-V sweeps, lock-in beta, fake data |
| `ivmodels` | single-level (Lorentzian) model at 0 K or finite temperature (closed digamma form, or the energy-grid Landauer integral of Rashid et al. 2025), Levenberg-Marquardt or bounded fits, most probable I-V curve, transition voltage spectroscopy |
| `events` | flickering and mechanical events in G(t) (Butterworth + Savitzky-Golay + derivative peaks), conductance-time histograms of holds |
| `egain` | large-area / EGaIn junctions: J from the contact area, shorts and yield, Gaussian log\|J\| per bias with junction-based confidence intervals, rectification, V_trans, beta and J0 from a length series |
| `figures`, `style`, `export` | plots, figure style, result files |

## Checks I did

- Conductance from the raw TDMS agrees with the G channel the acquisition software writes
  to within 0.003 decades (10^-4 to 10^-2 G0); piezo identical.
- Simulated data: peak at -3.54 (true -3.5), plateau length 0.42 nm (true 0.39 nm),
  n_TSE = 1.52 +/- 0.02 (true 1.5). With 5 % outlier windows Theil-Sen gives 1.40, OLS 1.33.
  dI/dV within 0.3 %, modulation beta 9.00 /nm (true 9.0).
- ADF test: same statistic, lag and p-value as statsmodels on 300 random series.
- Single-level model: on a model curve with eps0 = 0.8 eV and Gamma = 20 meV the fit gives 0.800 eV
  and 20.0 meV. TVS with both polarities gives eps0 = 0.800 eV and the asymmetry a = 0.100 (true 0.1).
- Junction yield on simulated data: 68.5 +/- 3.3 % (true 70 %) and 29.0 +/- 3.2 % (true 30 %).
- Finite temperature: closed form and energy-grid integral agree to 1e-9. With the Rashid et al. 2025 settings
  (300 K, grid, Levenberg-Marquardt) a model curve with eps0 = 1.10 eV, Gamma = 11.3 meV comes back exactly.
- Dwell time from the peak area: 135-146 ms (expected 140 ms on simulated data).
- Rashid et al. 2025 noise preset: n = 1.53 averaged over five simulated sets (true 1.5).
- EGaIn: beta = 0.904 +/- 0.016 (true 0.9), log R = 0.30 (true 0.3); over 20 repeats 65 % of the beta values
  fall within their 1-sigma error.

## Methods

- Flicker noise, bell: Adak et al., Nano Lett. 15, 4143 (2015), doi:10.1021/acs.nanolett.5b01270
- Flicker noise, zero correlation (OLS): Magyarkuti et al., Nanoscale 10, 3362 (2018), doi:10.1039/C7NR08354H
- Flicker noise, ADF + Theil-Sen: Morris et al., J. Phys. Chem. C 129, 4097 (2025), doi:10.1021/acs.jpcc.4c07780
- 2D correlation: Makk et al., ACS Nano 6, 3411 (2012), doi:10.1021/nn300440f
- Junction yield: Kamenetska et al., Phys. Rev. Lett. 102, 126803 (2009), doi:10.1103/PhysRevLett.102.126803
- Finite-temperature single-level fit on an energy grid, most probable I-V, noise protocol, events, dwell time:
  Rashid et al., J. Am. Chem. Soc. 147, 830 (2025), doi:10.1021/jacs.4c13531
- EGaIn statistics: Reus et al., J. Phys. Chem. C 116, 6714 (2012), doi:10.1021/jp210445y;
  beta and J0: Simeone et al., J. Am. Chem. Soc. 135, 18131 (2013), doi:10.1021/ja408652h
- Single-level model fit: Zotti et al., Small 6, 1529 (2010), doi:10.1002/smll.200902227
- Transition voltage spectroscopy: Beebe et al., Phys. Rev. Lett. 97, 026801 (2006), doi:10.1103/PhysRevLett.97.026801;
  Baldea, Phys. Rev. B 85, 035442 (2012), doi:10.1103/PhysRevB.85.035442
- Clustering: Lemmer et al., Nat. Commun. 7, 12922 (2016), doi:10.1038/ncomms12922;
  Cabosart et al., Appl. Phys. Lett. 114, 143102 (2019), doi:10.1063/1.5089198

Full reference list in the manual, chapter 15 (what to cite for which analysis: 15.7).

## Related open-source code

These were read while adding the EGaIn tab and the I-V and noise methods:

- [XMe_DataAnalysis](https://github.com/Pilab-XMU/XMe_DataAnalysis), Wenjing Hong's group (Xiamen University),
  Apache 2.0: break-junction, I-V, flicker-noise (PSD), correlation, spectral clustering, thermopower,
  electrochemistry and EGaIn modules. Please cite it as its README asks if you compare with it.
- [GaussFit](https://github.com/rchiechi/gaussfit), Chiechi group, GPL-3.0: EGaIn / CP-AFM statistics.
  Only its methods were used; no code was copied.

## License

Apache License 2.0, see [LICENSE](LICENSE). You can use, change and share the program, also
commercially; keep the license and copyright notices (see [NOTICE](NOTICE)). If you use it for a
publication, please cite the methods listed above and in chapter 15 of the manual.
