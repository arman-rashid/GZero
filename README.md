# GZero

**Break-junction analysis software.** GZero (G-zero, after the conductance quantum
G0 = 2e²/h) analyses single-molecule break-junction data (MCBJ and STM-BJ): raw TDMS, text
(.cvr/.txt/.csv) or Igor (.pxp/.ibw) files to conductance, trace detection, histograms, plateau
length, 2D correlation, clustering, flicker noise, I-V and piezo modulation.

The manual is `docs/User_Manual.pdf` (Help -> User manual in the program).

## Getting it

**Windows, no Python needed:** download `GZero-v1.1-windows-x64.zip` from the
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
| `origin` | rebuilding the plots in Origin (OriginLab) from the drawn matplotlib panels |
| `clustering` | 2D-histogram features, PCA, k-means / Gaussian mixture / Ward / spectral |
| `noise` | noise power, window selection, scaling exponent n (Theil-Sen, OLS, 2D Gaussian) |
| `stationarity` | ADF test, same results as statsmodels' adfuller |
| `iv`, `modulation`, `simulation` | I-V sweeps, lock-in beta, fake data |
| `figures`, `style`, `export` | plots, figure style, result files |

## Checks I did

- Conductance from the raw TDMS agrees with the G channel the acquisition software writes
  to within 0.003 decades (10^-4 to 10^-2 G0); piezo identical.
- Simulated data: peak at -3.54 (true -3.5), plateau length 0.42 nm (true 0.39 nm),
  n_TSE = 1.52 +/- 0.02 (true 1.5). With 5 % outlier windows Theil-Sen gives 1.40, OLS 1.33.
  dI/dV within 0.3 %, modulation beta 9.00 /nm (true 9.0).
- ADF test: same statistic, lag and p-value as statsmodels on 300 random series.

## Methods

- Flicker noise, bell: Adak et al., Nano Lett. 15, 4143 (2015), doi:10.1021/acs.nanolett.5b01270
- Flicker noise, zero correlation (OLS): Magyarkuti et al., Nanoscale 10, 3362 (2018), doi:10.1039/C7NR08354H
- Flicker noise, ADF + Theil-Sen: Morris et al., J. Phys. Chem. C 129, 4097 (2025), doi:10.1021/acs.jpcc.4c07780
- 2D correlation: Makk et al., ACS Nano 6, 3411 (2012), doi:10.1021/nn300440f
- Clustering: Lemmer et al., Nat. Commun. 7, 12922 (2016), doi:10.1038/ncomms12922;
  Cabosart et al., Appl. Phys. Lett. 114, 143102 (2019), doi:10.1063/1.5089198

Full reference list in the manual, chapter 15.

## License

Apache License 2.0, see [LICENSE](LICENSE). You can use, change and share the program, also
commercially; keep the license and copyright notices (see [NOTICE](NOTICE)). If you use it for a
publication, please cite the methods listed above and in chapter 15 of the manual.
