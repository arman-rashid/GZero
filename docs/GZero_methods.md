# GZero 1.3 to 1.5: new analysis methods

October 2026, arman-rashid

This note describes what changed in GZero 1.3, 1.4 and 1.5, why, and how it was checked. Version 1.3 adds three things that most break-junction papers report and GZero couldn't do yet: the junction yield with per-trace plateau statistics, single-level model fits of I-V curves, and transition voltage spectroscopy. Version 1.4 brings in the analysis from our own paper (Rashid et al., JACS 2025): the finite-temperature energy-grid fit, the noise protocol, event detection and the dwell time. It also adds two tabs for other kinds of measurements, EGaIn ensemble junctions and thermopower, and lets GZero read TDMS files of several GB. Version 1.5 adds change-point detection in breaking traces, the analysis of fixed-gap blinking (I(t)) measurements, six more file formats, and fixes for the problems found in a full review of the code.

## What's new

Each addition answers a question the existing tabs couldn't.

| Addition | Where | Question it answers | Module |
| --- | --- | --- | --- |
| Junction statistics | New tab after Plateau length | How many traces actually had a molecule, at what conductance, and did that drift during the run? | `gzero/junctions.py` |
| Single-level model fit | I-V tab | Where is the molecular level (eps0) and how strongly is it coupled to the electrodes (Gamma)? | `gzero/ivmodels.py` |
| Transition voltage spectroscopy | I-V tab | Where does the Fowler-Nordheim plot turn over, and what level position does that imply? | `gzero/ivmodels.py` |
| Events | New tab | Which traces flicker or switch mechanically? | `gzero/events.py` |
| EGaIn (ensemble) | New tab | What are the statistics of large-area junctions: yield, rectification and beta? | `gzero/egain.py` |
| Thermopower | New tab | What is the Seebeck coefficient, and is transport through the HOMO or the LUMO? | `gzero/thermo.py` |
| Large TDMS files | Data & Traces tab | Can a recording of more than 1 GB be analysed without running out of memory? | `gzero/streaming.py` |
| Change points (1.5) | New tab after Junction statistics | Where does each junction rupture, what is its own plateau conductance, and which traces have more than one step? | `gzero/changepoints.py` |
| Blinking (1.5) | New tab after Events | In a fixed-gap I(t) measurement: at what conductance do molecules bridge the gap, for how long, and how often? | `gzero/blinking.py` |
| More file formats (1.5) | Data & Traces, EGaIn, Thermopower | LabVIEW .lvm, Excel, MATLAB (also v7.3), NumPy, HDF5 and gzip-compressed text | `gzero/bj_io.py` |

One existing function was also fixed on the way. The averaged I-V curve used to take the overlap of all sweeps, so the half sweeps at the start and end of a recording could shrink it to almost nothing. It now only averages sweeps that cover the typical bias range.

## Papers reviewed

The methods sections of the most cited break-junction papers were searched for analyses that are routinely reported but that GZero didn't have. Four papers shaped version 1.3:

| Paper | What they do | What GZero does with it |
| --- | --- | --- |
| [Kamenetska et al., PRL 102, 126803 (2009)](https://doi.org/10.1103/PhysRevLett.102.126803) | Count how many breaking traces show a molecular plateau (junction formation probability) and how long the plateaus are, for a series of molecules | Junction yield per trace with a binomial error, plus yield and conductance in blocks of traces over the run |
| [Zotti et al., Small 6, 1529 (2010)](https://doi.org/10.1002/smll.200902227) | Fit MCBJ I-V curves of tolanes with thiol, nitro and cyano anchors with a single-level model to get eps0 and Gamma | Single-level fit of every curve and of the mean forward and backward curves, with errors and an optional bias asymmetry |
| [Beebe et al., PRL 97, 026801 (2006)](https://doi.org/10.1103/PhysRevLett.97.026801) | Plot ln(I/V^2) against 1/V and read off the minimum, the transition voltage V_t; it scales with the HOMO offset from photoemission | V_t for both polarities of every curve, with a Fowler-Nordheim plot of the mean curves |
| [Baldea, PRB 85, 035442 (2012)](https://doi.org/10.1103/PhysRevB.85.035442) | Shows that both polarities of V_t together give the level position eps0 and the bias asymmetry in closed form | eps0 and the asymmetry a from V_t+ and V_t-, so TVS and the fit can be compared curve by curve |

The 2D correlation, clustering and flicker-noise methods from earlier papers were already in GZero. Thermopower came later, in 1.4; shot noise is still out, because it needs hardware most setups don't have.

## Junction statistics

A trace counts as a molecular junction if it stays flat inside the molecular window for at least 0.1 nm, and the yield is the fraction of traces that do.

The tricky word is "flat". A trace that only tunnels still passes through a two-decade window in about 0.1-0.2 nm, so a plain length-in-window test would count almost every trace. GZero therefore looks at the local slope at each point:

```latex
s(z) = \frac{|\log G(z + w/2) - \log G(z - w/2)|}{w}
```

with w = 0.05 nm by default. Tunnelling falls by beta / ln 10 decades per nm, which is about 9 for gold in vacuum and 4-5 in solvents, while molecular plateaus usually fall by well under 1 decade per nm. A point counts as flat below 2 decades/nm, in between the two.

For every trace this gives a flat length (the summed distance of the flat points inside the window), a plateau conductance (the median log G of those points) and a yes or no: is the flat length at least *min plateau*? From these the tab works out the yield p with its binomial error sqrt(p(1-p)/N), the median plateau conductance and length, and the drift in decades per 1000 traces, which is a straight line through the medians of blocks of 100 traces. It plots the plateau conductance histogram, conductance and yield against trace number, the flat-length histogram with the threshold marked, and length against conductance.

In practice, Use peak from histogram fit sets the window to the fitted peak +/- 3 sigma. Then have a look at the flat-length histogram: tunnelling traces pile up near zero, molecular ones form a second hump, and the threshold belongs in the dip between them.

## Single-level model

The fit describes each I-V curve with one molecular level at energy eps0 from the Fermi level, broadened by its coupling Gamma to each electrode. At zero temperature the Landauer current through the resulting Lorentzian transmission has a closed form:

```latex
I(V) = N G_0 \Gamma \left[\arctan\frac{V/2 - \varepsilon(V)}{\Gamma} + \arctan\frac{V/2 + \varepsilon(V)}{\Gamma}\right], \quad \varepsilon(V) = \varepsilon_0 + aV
```

Energies are in eV and V in volts, so G0 times an energy in eV is directly a current in amperes. Both electrodes are coupled equally, and N, the number of molecules in parallel, stays at 1. The parameter a lets the level move with the bias in an asymmetric junction; it is 0 unless you tick *fit asymmetry a*. At low bias the model reduces to G/G0 = Gamma^2 / (eps0^2 + Gamma^2).

The current is first scaled by its largest value, so that every curve weighs the same. A least-squares fit is then started from four guesses of eps0 (0.4, 0.8, 1.5 and 3 times the largest bias), each with a Gamma that matches the measured low-bias conductance, and the one with the smallest residual wins. Errors come from the fit covariance, and R^2 is reported as well. This runs for every curve (up to 500), for the mean forward and backward curves, and for the most probable curve.

The model can't tell HOMO from LUMO transport, because both give the same I-V, so eps0 is reported as a positive number. The thermopower tab is the way to settle the sign.

## Transition voltage spectroscopy

The transition voltage V_t is where V^2/|I| peaks, which is the same as the minimum of the Fowler-Nordheim plot, ln(|I|/V^2) against 1/V. It needs no fit, and that's why so many papers report it.

For each curve and polarity GZero averages the current in 60 bias bins above 0.05 V, finds the largest V^2/|I| and refines it with a parabola through that point and its two neighbours. If the maximum sits at the end of the sweep, the bias never reached V_t, and the value is left empty rather than guessed.

In the single-level model with Gamma much smaller than eps0, the two polarities give the level position and the asymmetry exactly (Baldea 2012):

```latex
\varepsilon_0 = \frac{2|V_t^+ V_t^-|}{\sqrt{(V_t^+)^2 + (V_t^-)^2 + \frac{10}{3}|V_t^+ V_t^-|}}, \quad a = \frac{(V_t^+ + V_t^-)\,\varepsilon_0}{4|V_t^+ V_t^-|}
```

For a symmetric junction this reduces to eps0 = (sqrt(3)/2) V_t, about 0.87 V_t. The formula for a follows from the same stationarity condition, and it was checked numerically against the exact model (next section). Comparing the eps0 from TVS with the one from the fit is a quick way to see whether a single level describes the curves at all.

## Validation of the 1.3 methods

On data where the answer is known, every new method gets the true value back within its error. The SLM and TVS rows use model curves with 2 % noise; the yield rows use GZero's own simulator, going through the normal trace detection.

| Test | True value | GZero |
| --- | --- | --- |
| SLM fit, symmetric | eps0 0.800 eV, Gamma 10 meV | 0.801 eV, 10.0 meV (R^2 0.9997) |
| SLM fit, asymmetric | eps0 0.800 eV, a 0.100 | 0.795 eV, a 0.096 |
| TVS, symmetric | eps0 0.800 eV | V_t = +/-0.924 V, eps0 0.800 eV |
| TVS, asymmetric | eps0 0.800 eV, a 0.100 | V_t+ 1.191 V, V_t- -0.746 V, eps0 0.800 eV, a 0.100 |
| TVS, asymmetric | eps0 1.000 eV, a -0.080 | eps0 1.002 eV, a -0.080 |
| Self-test I-V | eps0 0.800 eV, Gamma 20 meV | fit 0.800 eV, 20.0 meV; TVS 0.801 eV |
| Junction yield, 200 traces | 70 % | 68.5 +/- 3.3 % |
| Junction yield, 200 traces | 30 % | 29.0 +/- 3.2 % |
| Plateau conductance | log G -3.5 | median -3.55 |

The small shift in the plateau conductance comes from the slope of the simulated plateaus (0.3 decades/nm), which pulls the median a little below where the plateau starts.

## Version 1.4: methods from Rashid et al. 2025

Every analysis in the SI of our paper, [Rashid et al., JACS 147, 830 (2025)](https://doi.org/10.1021/jacs.4c13531), is now in GZero. Only the DFT and transport calculations are not.

| SI section | What the paper does | In GZero |
| --- | --- | --- |
| 2.3 Flicker noise | Hold 160 ms at 100 mV, drop 10 ms, keep holds whose first and last 100 points sit within +/- 1 sigma of the peak, DFT squared, 100 Hz-1 kHz, n where Pearson r is smallest (0.3-2.3, step 0.01) | Flicker noise tab, preset "Rashid et al. 2025"; new *end points* setting; minimum-\|r\| grid value reported |
| 2.4 I-V fit | Landauer integral with Fermi functions on an energy grid (-10 to 10 eV, 0.02 meV), Gamma_L = Gamma_R, Levenberg-Marquardt fit of the most probable I-V | I-V tab, preset "Rashid et al. 2025"; energy grid or exact digamma form; most probable curve from a Gaussian per bias bin |
| 2.5 Mechanical modulation | 2D conductance-time histograms of holds with the most probable G per time | Events tab with segments = holds |
| 2.6 Event detection | Butterworth 200 Hz, Savitzky-Golay (order 20, 20 side points), derivative, peak thresholds 0.0007 / 0.002, width 50 | Events tab |
| 2.7 Dwell time | tau = A/w x 1000/f_s ms from the peak area | Histograms tab, with every peak fit |

At finite temperature the Landauer integral is done either on the energy grid we used in the paper or with the exact closed form, which uses the digamma function psi:

```latex
\int T(E)\, f(E - \mu)\, dE = \pi\Gamma \left[\frac{1}{2} - \frac{1}{\pi}\, \mathrm{Im}\, \psi\!\left(\frac{1}{2} + \frac{\Gamma - i(\mu - \varepsilon)}{2\pi k_B T}\right)\right]
```

The two agree to 10^-9. So the fits start from the closed form and the grid only polishes the result, which brings a grid fit down from 50 s to under 1 s.

| Test | True | GZero 1.4 |
| --- | --- | --- |
| Preset fit, model curve with the G1 values of Table S1, 300 K | eps0 1.10 eV, Gamma 11.3 meV | 1.100 eV, 11.30 meV |
| Same, 2 % noise, through the GUI (mean, most probable, single curves) | 1.10 eV, 11.3 meV | 1.098-1.100 eV, 11.3 meV |
| Noise preset, 5 simulated sets, holds on the plateau | n = 1.5 | 1.53 (Morris preset: 1.59) |
| Dwell time, 70 % molecular, 0.4 nm at 2 nm/s | 140 ms | 135-146 ms |

The SI doesn't give the temperature, so the preset assumes 300 K. Change it if you measured at another temperature.

## Events tab

The tab finds flickering and mechanical switching in G(t) the way SI 2.6 describes, and sorts every trace into mechanical, flickering or quiet. G(t), in G0, is smoothed with a 200 Hz, order 2 Butterworth low-pass, run forwards and backwards so the events don't shift in time, and then with a Savitzky-Golay filter of order 20 over 41 points. scipy's own coefficients collapse to about 10^-17 at that order, so GZero builds the filter from a QR decomposition in the Legendre basis instead, which matches scipy to 10^-14 at low order. The smoothed trace is differentiated with a central difference, and peaks of dG/dt above 0.0007 are flickers, above 0.002 mechanical events. Events while the junction is outside the molecular log G range are dropped as false positives.

The plots show an example trace with its events, its derivative with both thresholds, the share of each kind of trace, flickers per trace, and a before/after transition map. With holds you also get the conductance-time histogram with the most probable G.

Three details aren't pinned down in the SI, so GZero does what LabVIEW does by default. They're worth checking against the original VI:

- The 20 ms is LabVIEW's dt input, so dG/dt = (G[i+1] - G[i-1]) / (2 x 0.02 s). Since the dG/dt peak is about 2.2 times the step, the thresholds then correspond to steps of about 3 x 10^-4 G0 (flicker) and 9 x 10^-4 G0 (mechanical).
- "Width 50" is the window of the Peak Detector's quadratic fit, not a minimum peak width. Read as a minimum width at half height it would reject almost every step, because after the 200 Hz filter the derivative peaks are only 20-30 samples wide at 10 kHz. Both readings can be selected.
- Only rises are counted, as in the Peak Detector, so a flicker up and back counts once. *Count = both* counts every step.

On synthetic steps, a 4 x 10^-4 G0 flicker at 10^-3 G0 is found as a flicker, a jump from 10^-6 to 1.6 x 10^-3 G0 as mechanical, and a flat trace with 3 % noise gives nothing.

## EGaIn ensemble tab

Large-area junctions are analysed per junction, not per sweep: 20 junctions with 20 sweeps each are 20 independent measurements, not 400. The tab follows [Reus et al. 2012](https://doi.org/10.1021/jp210445y), as XMe and GaussFit do.

A file is either one junction, with the bias and current columns chosen in the settings, or many junctions stored as column pairs V, I, V, I, ..., which may differ in length. The current becomes J = I / A, with A from the contact diameter (as in XMe) or typed in. A junction whose \|J\| reaches the short threshold (100 A/cm2) is a short, and the yield is working junctions over all junctions. The sweeps go onto a common bias grid, and at each bias log10\|J\| is histogrammed and fitted with a Gaussian, which gives the centre and sigma_log. The confidence interval uses the number of junctions n_j as degrees of freedom, as GaussFit does. The tab also works out log R = log\|J(+V)\| - log\|J(-V)\| and V_trans per cycle, and fits the single-level model to the mean J(V) with N free. Across data sets with lengths d, it fits beta, weighting each point by the standard error of its mean ([Simeone et al. 2013](https://doi.org/10.1021/ja408652h) report 0.92 per carbon):

```latex
\mathrm{CI} = t_{1-\alpha/2,\; n_j - 1}\, \frac{\sigma_\mathrm{log}}{\sqrt{n_j - 1}}, \qquad \log_{10} |J| = \log_{10} J_0 - \frac{\beta\, d}{\ln 10}
```

| Check (simulated junctions) | True | Found |
| --- | --- | --- |
| beta, 5 sets of 30 junctions | 0.9 | 0.904 +/- 0.016 |
| Rectification log R at 0.5 V | 0.30 | 0.30 +/- 0.08 |
| Mean log\|J\| per bias | as simulated | inside its 95 % CI |
| beta error coverage, 20 repeats | 68 % within 1 sigma | 65 % |

The first version weighted beta by sigma_log and claimed errors three times too small. Weighting by the standard error fixed that, which is what the coverage row shows.

## Thermopower tab

The tab takes the thermoelectric voltage and the temperature (or temperature difference) of every sample and follows [Reddy et al. 2007](https://doi.org/10.1126/science.1137149) and [Widawsky et al. 2012](https://doi.org/10.1021/nl203634m). The temperature column can hold dT, or the hot-side temperature (then dT = T - T_ref), or be missing, in which case each data set gets its own dT. The samples, or the file means (one per junction), are grouped by dT, and a Gaussian is fitted to the thermovoltage histogram at every dT. A straight line through the most probable voltages gives S from its slope:

```latex
S = -\frac{dV}{d(\Delta T)} + S_\mathrm{electrode}
```

S > 0 means hole transport through the HOMO and S < 0 electron transport through the LUMO. S_i = -V_i/dT_i + S_electrode of every sample is histogrammed too, and plotted against log G when there is a conductance column. The error of each peak comes from a bootstrap over the files, because the samples within one trace aren't independent.

| Check (simulated junctions) | True | Found |
| --- | --- | --- |
| S, 160 traces at dT = 0-30 K | 12 uV/K | 11.9 +/- 0.4 uV/K |
| S, negative | -8 uV/K | -8.1 +/- 0.4 uV/K |
| Real files: V in mV, hot-side T, G window | -6 uV/K | -5.9 +/- 0.5 uV/K |
| Error coverage, 30 repeats | 68 % within 1 sigma | 80 % (slightly conservative) |

## Large TDMS files

Converting a 1 GB TDMS file in one go needs several GB of memory, since every step makes another array as long as the file. Files above 300 MB are therefore read in pieces of 2 minutes. Loading keeps only a min/max overview for the plot. Cutting the traces reads each piece with a 60 s margin on both sides and keeps only the traces and holds that start inside the piece, so a trace crossing a border is found once and in full. On a simulated recording this gives exactly the same traces as reading the whole file, and the self-test checks that. A 1.29 GB file (90 minutes at 10 kHz) now loads in 2.8 s and needs 0.7 GB of memory at most, against 15.6 s and 4.4 GB before.

## Version 1.5: change points, blinking, file formats

For 1.5 the recent literature on break-junction analysis was searched for methods that answer questions the existing tabs can't. Two stood out. Change-point detection gives one conductance per junction, measured over its own plateau, and finds traces with more than one step. Blinking (the I(t) method) is a whole measurement mode in its own right, and GZero had no tab for it.

### Change points

The rupture of each trace follows [Hamill et al. 2024](https://arxiv.org/abs/2401.12769). Inside the window between log G = -2 and -5.5, log G against distance is fitted by one straight line and by two lines meeting at a break k. The Chow statistic ([Chow 1960](https://doi.org/10.2307/1910133)) compares the two:

```latex
F(k) = \frac{\left(S_1 - S_2(k)\right)/p}{S_2(k)/(n - 2p)}, \qquad p = 2
```

The k with the largest F is the rupture, and the median log G before it is the plateau conductance of that junction. A large F alone isn't enough, because a plain tunnelling decay also bends slightly and its noise is small. So a trace only counts as a junction if the line before the break is also flat: |slope| below 2 decades/nm, the same criterion as the Junction statistics tab. On simulated traces the slope before the break is -0.3 decades/nm for every molecular trace and -4.3 for every blank one, so the two groups don't overlap.

The steps use PELT ([Killick et al. 2012](https://doi.org/10.1080/01621459.2012.737745)) for a piecewise-constant fit of log G with a penalty of 2 ln(n) per change point and the noise sigma estimated from the trace (at least 0.05 decades). A variable number of change points per trace is the idea of PELT-R ([Vaidya, Dyer and Monti 2026](https://meetings-archive.aps.org/smt/2026/mar-s71/8/)). A tunnelling ramp is cut into many short levels by such a fit, so only levels that are themselves flat and at least 0.05 nm long count as plateau levels.

| Check | True | Found |
| --- | --- | --- |
| Rupture, 200 traces, 70 % molecular | 70 % | 70.0 % junctions |
| Rupture, 200 traces, 30 % molecular | 30 % | 31.5 % |
| Plateau conductance per junction | log G -3.5 (spread 0.25) | median -3.55, spread 0.24 |
| PELT against an exhaustive search, 5 random series | same change points | identical |
| PELT, trace with plateaus at -2.9 and -3.9 (0.3 nm each) | two levels | -2.90 and -3.90, 0.30 and 0.32 nm |

### Blinking

With the gap held fixed, a molecule bridging it shows up as a jump of the current to the molecular level, and its leaving as a jump back ([Haiss et al. 2003](https://doi.org/10.1021/ja038214e), [2006](https://doi.org/10.1038/nmat1781)). The baseline is a low percentile (20 %) of log G in 0.5 s blocks, so it follows drift without being pulled up by the blinks. A blink starts 0.6 decades above the baseline and ends below 0.3 decades (hysteresis). The blink conductance is a Gaussian fit to the blink levels. The on and off times are fitted with the maximum-likelihood lifetime of an exponential above the shortest time measured:

```latex
\tau = \frac{N}{N - 1}\,\left(\bar{t} - t_\mathrm{min}\right), \qquad \delta\tau = \tau / \sqrt{N}
```

| Check (simulated 20 s recordings) | True | Found |
| --- | --- | --- |
| Blinks found | 143, 172, 186 | 137, 171, 179 |
| Blink conductance | log G -3.5 | -3.51 to -3.53 (+/- 0.01) |
| tau_on | 20 ms | 19.7-23.4 ms (+/- 1.5-2.0) |
| tau_off | 100 ms | 91-121 ms (+/- 7-10) |

### File formats

Every non-TDMS format is turned into the same list of named columns, so the column settings work unchanged. New in 1.5: LabVIEW `.lvm` (header names, `Delta_X` as the sampling interval), Excel `.xlsx` (every sheet), MATLAB `.mat` (struct fields, v7.3 through h5py), NumPy `.npy` / `.npz`, HDF5 `.h5` / `.hdf5` (data sets by path, sampling rate from an attribute), and gzip-compressed text. The self-test writes a model trace in each of these formats, reads it back and checks that the conductance agrees to 10^-5 decades.

### Fixes

A complete review of the code found these problems, all fixed in 1.5:

| Problem | Effect before 1.5 |
| --- | --- |
| The delimiter of a text file was guessed from the first line, even if that was a header with a comma | Such files could not be read |
| Missing current samples (NaN) and zero bias became G = 100 G0 | Fake metal-contact spikes; now they go to the floor |
| The diode calibration left NaN inputs unset | Random values in the converted current |
| The noise power of each window was integrated after interpolation onto the frequency grid of the first window | With whole-hold windows of different length, NP depended on which hold came first; now each window is integrated on its own frequencies |
| A peak fit was kept when the histogram settings changed but the number of traces didn't | An old fit and dwell time could be shown and saved with a new histogram |
| Junction statistics, Correlation, Flicker noise, Events, I-V and Piezo modulation kept their results after Clear or new traces | Old results could be saved under the name of new data |
| Plateau export with no valid window | Error during Save results |
| The largest point fell out of the binned noise trend | Last bin one point short |
| The 2D histogram in `density` mode wasn't divided by the bin area | Wrong scale |
| A recording with fewer than 3 piezo samples | Error during trace detection |
| The bell fit was repeated at every n of the scan (up to 200 fits) | Very slow with few windows; now at most 40 fits, interpolated |


## Open-source code reviewed

Two public codes cover the same ground, and both were read in full for their EGaIn and I-V methods.

| Code | Group | License | Modules | How GZero used it |
| --- | --- | --- | --- | --- |
| [XMe_DataAnalysis](https://github.com/Pilab-XMU/XMe_DataAnalysis) | Wenjing Hong, Xiamen University | Apache 2.0 | Basic histograms, correlation, spectral clustering, PSD noise, I-V, single trace, thermopower, electrochemistry, EGaIn | EGaIn: J from contact diameter, log\|J\| per bias with Gaussian fits. The methods were re-implemented; XMe asks to be cited when compared with |
| [GaussFit](https://github.com/rchiechi/gaussfit) | Ryan Chiechi | GPL-3.0 | EGaIn / CP-AFM parsing, Gaussian log\|J\|, rectification, V_trans, dJ/dV, lag filter, SLM | Junction-based CI, rectification and V_trans per cycle. Methods only: GPL code can't go into Apache-licensed GZero |

Of the other Hong-group analyses in XMe, GZero already had correlation, spectral clustering and flicker noise, and thermopower is in now. Electrochemical gating needs a measurement type GZero doesn't read yet.

## Sources

Code: [arman-rashid/GZero](https://github.com/arman-rashid/GZero), release v1.5. The full method description is in the user manual (`docs/User_Manual.pdf`), sections 3.1, 4.7, 5.2, 6.4, 6.5, 9.5, 9.8, 9.9 and 10.1-10.6, references 26-43.

- M. Kamenetska et al., "Formation and evolution of single-molecule junctions", Phys. Rev. Lett. 102, 126803 (2009). [doi:10.1103/PhysRevLett.102.126803](https://doi.org/10.1103/PhysRevLett.102.126803)
- L. A. Zotti et al., "Revealing the role of anchoring groups in the electrical conduction through single-molecule junctions", Small 6, 1529-1535 (2010). [doi:10.1002/smll.200902227](https://doi.org/10.1002/smll.200902227)
- J. M. Beebe et al., "Transition from direct tunneling to field emission in metal-molecule-metal junctions", Phys. Rev. Lett. 97, 026801 (2006). [doi:10.1103/PhysRevLett.97.026801](https://doi.org/10.1103/PhysRevLett.97.026801)
- I. Baldea, "Ambipolar transition voltage spectroscopy: analytical results and experimental agreement", Phys. Rev. B 85, 035442 (2012). [doi:10.1103/PhysRevB.85.035442](https://doi.org/10.1103/PhysRevB.85.035442)
- U. Rashid et al., "Mapping the extended ground state reactivity landscape of a photoswitchable molecule at a single molecular level", J. Am. Chem. Soc. 147, 830-840 (2025). [doi:10.1021/jacs.4c13531](https://doi.org/10.1021/jacs.4c13531) (SI sections 2.3-2.7)
- L. J. Reus et al., "Statistical tools for analyzing measurements of charge transport", J. Phys. Chem. C 116, 6714 (2012). [doi:10.1021/jp210445y](https://doi.org/10.1021/jp210445y)
- F. C. Simeone et al., "Defining the value of injection current and effective electrical contact area for EGaIn-based molecular tunneling junctions", J. Am. Chem. Soc. 135, 18131 (2013). [doi:10.1021/ja408652h](https://doi.org/10.1021/ja408652h)
- P. Reddy et al., "Thermoelectricity in molecular junctions", Science 315, 1568 (2007). [doi:10.1126/science.1137149](https://doi.org/10.1126/science.1137149)
- J. R. Widawsky et al., "Simultaneous determination of conductance and thermopower of single molecule junctions", Nano Lett. 12, 354 (2012). [doi:10.1021/nl203634m](https://doi.org/10.1021/nl203634m)
- E. M. Opodi et al., "A single level tunneling model for molecular junctions: evaluating the simulation methods", Phys. Chem. Chem. Phys. 24, 11958 (2022). [doi:10.1039/D1CP05807J](https://doi.org/10.1039/D1CP05807J)
- J. M. Hamill et al., "Improving single-molecule conductance measurements with change point detection from the econometrics toolbox", arXiv:2401.12769 (2024). [arxiv.org/abs/2401.12769](https://arxiv.org/abs/2401.12769)
- G. C. Chow, "Tests of equality between sets of coefficients in two linear regressions", Econometrica 28, 591 (1960). [doi:10.2307/1910133](https://doi.org/10.2307/1910133)
- R. Killick, P. Fearnhead, I. A. Eckley, "Optimal detection of changepoints with a linear computational cost", J. Am. Stat. Assoc. 107, 1590 (2012). [doi:10.1080/01621459.2012.737745](https://doi.org/10.1080/01621459.2012.737745)
- A. Vaidya, J. Dyer, O. L. A. Monti, "Understanding conductance fluctuations in single-molecule break junctions: a statistical approach", APS March Meeting 2026 (abstract). [meetings-archive.aps.org](https://meetings-archive.aps.org/smt/2026/mar-s71/8/)
- W. Haiss et al., "Redox state dependence of single molecule conductivity", J. Am. Chem. Soc. 125, 15294 (2003). [doi:10.1021/ja038214e](https://doi.org/10.1021/ja038214e)
- W. Haiss et al., "Precision control of single-molecule electrical junctions", Nature Materials 5, 995 (2006). [doi:10.1038/nmat1781](https://doi.org/10.1038/nmat1781)
- [XMe_DataAnalysis](https://github.com/Pilab-XMU/XMe_DataAnalysis) and [GaussFit](https://github.com/rchiechi/gaussfit) on GitHub
