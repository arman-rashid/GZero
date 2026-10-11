# GZero 1.3 and 1.4: new analysis methods

October 2026, arman-rashid

GZero 1.3 adds three methods that are standard in break-junction papers: the junction yield with per-trace plateau statistics, single-level model fits of I-V curves, and transition voltage spectroscopy. GZero 1.4 adds the analysis methods of Rashid et al., JACS 2025 (finite-temperature energy-grid fit, noise protocol, event detection, dwell time) and an EGaIn ensemble tab built after the Hong (XMe) and Chiechi (GaussFit) codes.

## What's new

Each new method answers a question the existing tabs couldn't. All of them run on data GZero already loads, so no new file formats are needed.

| Addition | Where | Question it answers | Module |
| --- | --- | --- | --- |
| Junction statistics | New tab after Plateau length | How many traces actually had a molecule, at what conductance, and did that drift during the run? | `gzero/junctions.py` |
| Single-level model fit | I-V tab | Where is the molecular level (eps0) and how strongly is it coupled to the electrodes (Gamma)? | `gzero/ivmodels.py` |
| Transition voltage spectroscopy | I-V tab | Where does the Fowler-Nordheim plot turn over, and what level position does that imply? | `gzero/ivmodels.py` |
| Events | New tab | Which traces flicker or switch mechanically? | `gzero/events.py` |
| EGaIn (ensemble) | New tab | Statistics of large-area junctions, yield, rectification and beta | `gzero/egain.py` |
| Thermopower | New tab | What is the Seebeck coefficient, and is transport through the HOMO or the LUMO? | `gzero/thermo.py` |

One existing function got a fix along the way. The averaged I-V curve used to take the overlap of all sweeps, so the half sweeps at the start and end of a recording could shrink it to nothing. It now averages only sweeps that cover the typical bias range.

## Papers reviewed

I went through the methods sections of commonly cited break-junction papers and looked for analyses people report routinely that GZero didn't have yet. Four papers define what was added in 1.3:

| Paper | What they do | What GZero does with it |
| --- | --- | --- |
| [Kamenetska et al., PRL 102, 126803 (2009)](https://doi.org/10.1103/PhysRevLett.102.126803) | Count how many breaking traces show a molecular plateau (junction formation probability) and how long the plateaus are, for a series of molecules | Junction yield per trace with a binomial error, plus yield and conductance in blocks of traces over the run |
| [Zotti et al., Small 6, 1529 (2010)](https://doi.org/10.1002/smll.200902227) | Fit MCBJ I-V curves of tolanes with thiol, nitro and cyano anchors with a single-level model to get eps0 and Gamma | Single-level fit of every curve and of the mean forward and backward curves, with errors and an optional bias asymmetry |
| [Beebe et al., PRL 97, 026801 (2006)](https://doi.org/10.1103/PhysRevLett.97.026801) | Plot ln(I/V^2) against 1/V and read off the minimum, the transition voltage V_t; it scales with the HOMO offset from photoemission | V_t for both polarities of every curve, with a Fowler-Nordheim plot of the mean curves |
| [Baldea, PRB 85, 035442 (2012)](https://doi.org/10.1103/PhysRevB.85.035442) | Shows that both polarities of V_t together give the level position eps0 and the bias asymmetry in closed form | eps0 and the asymmetry a from V_t+ and V_t-, so TVS and the fit can be compared curve by curve |

The 2D correlation, clustering and flicker-noise methods from earlier papers were already in GZero. I left out methods that need hardware most setups don't have, like thermopower and shot noise.

## Junction statistics

A trace counts as a molecular junction if it stays flat inside the molecular window for at least 0.1 nm. The yield is the fraction of traces that do.

The tricky part is "flat". A trace that only tunnels still passes through a two-decade window in about 0.1-0.2 nm, so a plain length-in-window test counts almost every trace. GZero therefore checks the local slope at each point:

```latex
s(z) = \frac{|\log G(z + w/2) - \log G(z - w/2)|}{w}
```

w is 0.05 nm by default. Tunnelling falls at beta / ln 10 decades per nm: about 9 for gold in vacuum and 4-5 in solvents. Molecular plateaus usually fall by well under 1 decade per nm. A point counts as flat below 2 decades/nm, which sits between the two.

Per trace you get three numbers:

- **flat length**: summed distance of the flat points inside the window
- **plateau conductance**: median log G of those points
- **junction yes/no**: flat length of at least *min plateau*

From these the tab gives the yield p with its binomial error sqrt(p(1-p)/N), the median plateau conductance and length, and a drift in decades per 1000 traces. The drift is a straight line fitted to the block medians, with 100 traces per block by default. The five plots are the plateau conductance histogram, conductance against trace number, yield against trace number, the flat-length histogram with the threshold marked, and length against conductance.

Use peak from histogram fit sets the window to the fitted peak +/- 3 sigma. Then look at the flat-length histogram. Tunnelling traces pile up near zero and molecular ones form a second hump, and the threshold belongs in the dip between them.

## Single-level model

The fit describes each I-V curve with one molecular level at energy eps0 from the Fermi level, broadened by the coupling Gamma to each electrode. At zero temperature the Landauer current through a Lorentzian transmission has a closed form:

```latex
I(V) = N G_0 \Gamma \left[\arctan\frac{V/2 - \varepsilon(V)}{\Gamma} + \arctan\frac{V/2 + \varepsilon(V)}{\Gamma}\right], \quad \varepsilon(V) = \varepsilon_0 + aV
```

Energies are in eV and V in volts, so G0 times an energy in eV is directly a current in amperes. Both electrodes are coupled equally. N is the number of molecules in parallel and stays at 1. a lets the level move with the bias in an asymmetric junction; it is 0 unless you tick *fit asymmetry a*. At low bias the model gives G/G0 = Gamma^2 / (eps0^2 + Gamma^2).

How the fit works:

1. The current is scaled by its largest value so all curves weigh the same.
2. A least-squares fit starts from four guesses of eps0 (0.4, 0.8, 1.5 and 3 times the largest bias). Each guess gets a Gamma that matches the measured low-bias conductance.
3. The fit with the smallest residual is kept. Errors come from the fit covariance, and R^2 is reported.
4. This is done for every curve (up to 500), for the mean forward and backward curves, and for the most probable curve.

The model can't tell HOMO from LUMO transport, because both give the same I-V. eps0 is therefore reported as a positive number.

## Transition voltage spectroscopy

The transition voltage V_t is where V^2/|I| peaks, which is the same as the minimum of the Fowler-Nordheim plot, ln(|I|/V^2) against 1/V. It needs no fit, which is why so many papers report it.

For each curve and polarity GZero averages the current in 60 bias bins above 0.05 V and finds the largest V^2/|I|. It then refines that point with a parabola through it and its two neighbours. If the maximum sits at the end of the sweep, the bias never reached V_t, and the value is left empty rather than guessed.

In the single-level model with Gamma much smaller than eps0, the two polarities give the level position and the asymmetry exactly (Baldea 2012):

```latex
\varepsilon_0 = \frac{2|V_t^+ V_t^-|}{\sqrt{(V_t^+)^2 + (V_t^-)^2 + \frac{10}{3}|V_t^+ V_t^-|}}, \quad a = \frac{(V_t^+ + V_t^-)\,\varepsilon_0}{4|V_t^+ V_t^-|}
```

For a symmetric junction this reduces to eps0 = (sqrt(3)/2) V_t, about 0.87 V_t. The a formula follows from the same stationarity condition and was checked numerically against the exact model (next section). Comparing the eps0 from TVS with the eps0 from the fit is a quick test of whether a single level describes the curves at all.

## Validation of the 1.3 methods

On data with known answers, every new method recovers the true value within its error. The SLM and TVS rows are model curves with 2 % noise. The yield rows are GZero's own simulator going through the normal trace detection.

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

The slight shift of the plateau conductance comes from the simulated plateau slope (0.3 decades/nm), which pulls the median a little below the plateau's start.

## Version 1.4: methods from Rashid et al. 2025

Every analysis section of the SI of [Rashid et al., JACS 147, 830 (2025)](https://doi.org/10.1021/jacs.4c13531) is now in GZero, except the DFT and transport calculations.

| SI section | What the paper does | In GZero |
| --- | --- | --- |
| 2.3 Flicker noise | Hold 160 ms at 100 mV, drop 10 ms, keep holds whose first and last 100 points sit within +/- 1 sigma of the peak, DFT squared, 100 Hz-1 kHz, n where Pearson r is smallest (0.3-2.3, step 0.01) | Flicker noise tab, preset "Rashid et al. 2025"; new *end points* setting; minimum-\|r\| grid value reported |
| 2.4 I-V fit | Landauer integral with Fermi functions on an energy grid (-10 to 10 eV, 0.02 meV), Gamma_L = Gamma_R, Levenberg-Marquardt fit of the most probable I-V | I-V tab, preset "Rashid et al. 2025"; energy grid or exact digamma form; most probable curve from a Gaussian per bias bin |
| 2.5 Mechanical modulation | 2D conductance-time histograms of holds with the most probable G per time | Events tab with segments = holds |
| 2.6 Event detection | Butterworth 200 Hz, Savitzky-Golay (order 20, 20 side points), derivative, peak thresholds 0.0007 / 0.002, width 50 | Events tab |
| 2.7 Dwell time | tau = A/w x 1000/f_s ms from the peak area | Histograms tab, with every peak fit |

At finite temperature the Landauer integral is evaluated either on the energy grid of the paper or with the exact closed form, which uses the digamma function psi:

```latex
\int T(E)\, f(E - \mu)\, dE = \pi\Gamma \left[\frac{1}{2} - \frac{1}{\pi}\, \mathrm{Im}\, \psi\!\left(\frac{1}{2} + \frac{\Gamma - i(\mu - \varepsilon)}{2\pi k_B T}\right)\right]
```

The two agree to 10^-9. Fits therefore start from the closed form and the grid only polishes, which brings a grid fit down from 50 s to under 1 s.

| Test | True | GZero 1.4 |
| --- | --- | --- |
| Preset fit, model curve with the G1 values of Table S1, 300 K | eps0 1.10 eV, Gamma 11.3 meV | 1.100 eV, 11.30 meV |
| Same, 2 % noise, through the GUI (mean, most probable, single curves) | 1.10 eV, 11.3 meV | 1.098-1.100 eV, 11.3 meV |
| Noise preset, 5 simulated sets, holds on the plateau | n = 1.5 | 1.53 (Morris preset: 1.59) |
| Dwell time, 70 % molecular, 0.4 nm at 2 nm/s | 140 ms | 135-146 ms |

The SI gives no temperature, so the preset assumes 300 K. Change it if your measurements ran at another temperature.

## Events tab

The tab finds flickering and mechanical switching in G(t) with the four steps of SI 2.6, and classifies each trace as mechanical, flickering or quiet.

1. Butterworth low-pass, 200 Hz, order 2, applied forwards and backwards so events aren't shifted in time.
2. Savitzky-Golay smoothing, order 20 over 41 points. scipy's own coefficients collapse to about 10^-17 at this order, so GZero builds the filter from a Legendre-basis QR. It matches scipy to 10^-14 at low order.
3. Second-order central derivative dG/dt with G in G0.
4. Peaks above 0.0007 are flickers and above 0.002 mechanical events. Events outside the molecular log G range are dropped as false positives.

The plots are an example trace with the events marked, its derivative with both thresholds, the share of each trace class, flickers per trace, and a before/after transition map. With holds you also get the conductance-time histogram with the most probable G.

**Three choices to confirm.** The SI doesn't pin these down, so GZero copies LabVIEW's defaults:

- [ ] **dt = 20 ms is LabVIEW's dt input**, so dG/dt = (G[i+1] - G[i-1]) / (2 x 0.02 s). The thresholds then correspond to steps of about 3 x 10^-4 G0 (flicker) and 9 x 10^-4 G0 (mechanical), since the dG/dt peak is about 2.2 times the step.
- [ ] **"Width 50" is the Peak Detector's quadratic-fit window**, not a minimum peak width. Read as a minimum width at half height, it rejects almost every step, because 200 Hz-filtered steps give derivative peaks only 20-30 samples wide at 10 kHz. Both readings are selectable.
- [ ] **Only rises are counted** (Peak Detector default), so a flicker up and back counts once. *Count = both* counts every step.

On synthetic steps a 4 x 10^-4 G0 flicker at 10^-3 G0 is found as a flicker, a jump from 10^-6 to 1.6 x 10^-3 G0 as mechanical, and a flat trace with 3 % noise gives nothing.

## EGaIn ensemble tab

Large-area junctions are analysed per junction, not per sweep. 20 junctions with 20 sweeps each are 20 independent measurements, not 400. The tab follows [Reus et al. 2012](https://doi.org/10.1021/jp210445y), as XMe and GaussFit also do.

1. Each file is one junction. J = I / A, with A from the contact diameter (XMe's approach) or typed in.
2. A junction whose \|J\| reaches the short threshold (100 A/cm2) is a short. The yield is working junctions over all junctions.
3. Sweeps go onto a common bias grid. At each bias, log10\|J\| is histogrammed and fitted with a Gaussian to give the centre and sigma_log.
4. The confidence interval uses the number of junctions n_j as degrees of freedom (GaussFit's convention); formula below.
5. Also computed: log R = log\|J(+V)\| - log\|J(-V)\| per cycle, V_trans per cycle, and a single-level fit of the mean J(V) with free N.
6. Fit beta across data sets with typed-in lengths d (formula below). Each point is weighted by the standard error of its mean ([Simeone et al. 2013](https://doi.org/10.1021/ja408652h) report 0.92 per carbon).

```latex
\mathrm{CI} = t_{1-\alpha/2,\; n_j - 1}\, \frac{\sigma_\mathrm{log}}{\sqrt{n_j - 1}}, \qquad \log_{10} |J| = \log_{10} J_0 - \frac{\beta\, d}{\ln 10}
```

| Check (simulated junctions) | True | Found |
| --- | --- | --- |
| beta, 5 sets of 30 junctions | 0.9 | 0.904 +/- 0.016 |
| Rectification log R at 0.5 V | 0.30 | 0.30 +/- 0.08 |
| Mean log\|J\| per bias | as simulated | inside its 95 % CI |
| beta error coverage, 20 repeats | 68 % within 1 sigma | 65 % |

The first version weighted beta by sigma_log and claimed errors 3 times too small. Weighting by the standard error fixed that, which is what the coverage row shows.

**File layouts.** A file is either one junction (bias and current columns chosen in the settings) or many junctions as column pairs V, I, V, I, ... Pairs may differ in length.

## Thermopower tab

The tab takes the thermoelectric voltage and the temperature (or temperature difference) of every sample and follows [Reddy et al. 2007](https://doi.org/10.1126/science.1137149) and [Widawsky et al. 2012](https://doi.org/10.1021/nl203634m):

1. Group the samples (or the file means, one per junction) by dT and fit a Gaussian to the thermovoltage histogram at every dT.
2. Fit the most probable voltage against dT with a straight line; the slope gives S (formula below). S > 0 means hole transport (HOMO), S < 0 electron transport (LUMO).
3. Histogram S_i = -V_i/dT_i + S_electrode of every sample, and plot it against log G when a conductance column is given.

```latex
S = -\frac{dV}{d(\Delta T)} + S_\mathrm{electrode}
```

The error of each peak comes from a bootstrap over the files, because the samples of one trace aren't independent. The temperature column can hold dT, the hot-side temperature (dT = T - T_ref), or be absent, in which case each data set gets its own dT.

| Check (simulated junctions) | True | Found |
| --- | --- | --- |
| S, 160 traces at dT = 0-30 K | 12 uV/K | 11.9 +/- 0.4 uV/K |
| S, negative | -8 uV/K | -8.1 +/- 0.4 uV/K |
| Real files: V in mV, hot-side T, G window | -6 uV/K | -5.9 +/- 0.5 uV/K |
| Error coverage, 30 repeats | 68 % within 1 sigma | 80 % (slightly conservative) |

## Open-source code reviewed

Two public codes cover the same ground and were read in full for their EGaIn and I-V methods.

| Code | Group | License | Modules | How GZero used it |
| --- | --- | --- | --- | --- |
| [XMe_DataAnalysis](https://github.com/Pilab-XMU/XMe_DataAnalysis) | Wenjing Hong, Xiamen University | Apache 2.0 | Basic histograms, correlation, spectral clustering, PSD noise, I-V, single trace, thermopower, electrochemistry, EGaIn | EGaIn: J from contact diameter, log\|J\| per bias with Gaussian fits. The methods were re-implemented; XMe asks to be cited when compared with |
| [GaussFit](https://github.com/rchiechi/gaussfit) | Ryan Chiechi | GPL-3.0 | EGaIn / CP-AFM parsing, Gaussian log\|J\|, rectification, V_trans, dJ/dV, lag filter, SLM | Junction-based CI, rectification and V_trans per cycle. Methods only: GPL code can't go into Apache-licensed GZero |

GZero already covers the other Hong-group analyses in XMe: correlation, spectral clustering and flicker noise. Thermopower and electrochemical gating need measurement types GZero doesn't read yet.

## Sources

Code: [arman-rashid/GZero](https://github.com/arman-rashid/GZero), release v1.4. The full method description is in the user manual (`docs/User_Manual.pdf`), sections 5.2, 6.4, 9.5, 9.8 and 10.2-10.6, references 26-37.

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
- [XMe_DataAnalysis](https://github.com/Pilab-XMU/XMe_DataAnalysis) and [GaussFit](https://github.com/rchiechi/gaussfit) on GitHub
