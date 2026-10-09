# GZero - User Manual

Break-junction analysis software. Version 1.1, October 2026

## 1. Introduction and installation

I wrote GZero to analyse single-molecule break-junction data (MCBJ and STM-BJ) in one place, from the raw TDMS files of the setup to the figures for a paper. It converts the raw signals to conductance, cuts out the single traces, makes the usual histograms, measures plateau lengths, calculates 2D correlation maps, clusters traces and gets the flicker-noise exponent. I-V sweeps and piezo modulation are covered too, and there is a simulator that makes fake data with known answers, which I use to check settings.

### 1.1 The tabs

| Tab | What it's for |
| --- | --- |
| Convert TDMS | Raw TDMS channels to conductance, `.cvr` export |
| Data & Traces | Loading files, finding and cutting traces, looking through them, marking bad ones |
| Histograms | 1D log G histogram, 2D conductance-distance histogram, Gaussian fits |
| Plateau length | How long the molecular plateau is, in up to three windows; tunnelling decay |
| Correlation | 2D cross-correlation of conductance values within traces |
| Clustering | PCA and clustering of traces into groups |
| Flicker noise | Noise power against conductance and the scaling exponent n |
| I-V | Forward and backward bias sweeps, dI/dV, 2D I-V histograms |
| Piezo modulation | Decay constant beta from a modulated piezo |
| Simulation | Synthetic recordings |

### 1.2 What you need

Windows 10 or 11 (64 bit) for the exe. From Python it should also run on macOS and Linux, though I've only used it on Windows. 8 GB of RAM is fine for normal files. Big files are loaded completely into memory: my 1 GB test file (73 min at 10 kHz) took about 7 s to load, and I'd plan roughly 3-4 GB of RAM per GB of TDMS. The window needs about 1400 x 900 pixels to be comfortable.

### 1.3 Installing and starting

The simplest way is the exe. Copy the whole `GZero` folder somewhere, e.g. `C:\Programs\GZero`, and double-click `GZero.exe`. Don't move the exe out of the folder, it needs the files next to it. This manual is in `docs\User_Manual.pdf` and also opens from Help -> User manual (F1).

The source code is on GitHub: https://github.com/arman-rashid/GZero (Apache License 2.0). To run it from Python (3.9 or newer), go to the project folder (the one with `GZero.py`) and run

```
py -m pip install -r requirements.txt
py GZero.py
```

It needs PySide6, numpy, scipy, matplotlib, pandas, scikit-learn, npTDMS and igor2 (for Igor `.pxp` / `.ibw` files).

If you change the code and want a new exe: `py -m pip install pyinstaller`, then `py build_exe.py` in the same folder. The new version ends up in `dist\GZero\`.

### 1.4 The window

Every tab has the settings and buttons on the left and the plots on the right. The toolbar above each plot lets you zoom, pan and save what's on screen. When you save results, each plot is also written as a separate figure in journal format (13.1). Most settings have a tooltip, so hover over them if something isn't clear. The status bar at the bottom says what the program is doing; longer calculations run in the background and the window doesn't freeze.

Settings are only kept while the program is open. File -> Save settings writes them to a file and File -> Load settings reads them back. You can also load the `*_config.json` from an old results folder, which is the easiest way to repeat an analysis with exactly the same settings.

## 2. Quick start

This is the short version of a normal analysis, from a raw file to a histogram and a noise exponent. The chapter numbers point to the details.

1. Go to Data & Traces, click Load files and pick a `.tdms`, `.cvr`, text or Igor `.pxp` file (several at once is fine, or use Load folder). For everything except TDMS check the column settings first: which column is what, and the bias if the file has no voltage column (4.2).
2. The traces are found automatically after the first load. Under the button you see how many were found, how many were accepted and why the rest were thrown out.
3. Click through a few traces with < and > or the arrow keys. Green means accepted, red rejected. If good-looking traces are being rejected, change the limits (4.4) and press Detect / cut traces again.
4. G and B mark the current trace as good or bad by hand.
5. Open Histograms. The 1D and 2D histograms are already there.
6. Put *fit from* and *fit to* around the molecular peak, choose how many Gaussians, click Fit peaks (5).
7. In Plateau length, set window 1 around the peak and click Calculate (6).
8. In Flicker noise, click the Morris et al. 2025 preset, then Peak from histogram fit, then Run noise analysis. n and its error appear under the button (9).
9. Correlation and Clustering if you need them (7, 8).
10. File -> Save results (Ctrl+S), choose where and give the folder a name. Everything goes in there: text files, figures, settings (13).

If you just want to try the program, go to Simulation and click Generate. You get a fake data set where you know the right answers, and you can see if the analysis finds them (12).

## 3. Data formats and the Convert tab

You only need the Convert tab for raw `.tdms` files. It turns them into conductance (G/G0) and either writes `.cvr` files or loads the result straight into the analysis. All other file types are loaded in Data & Traces (4).

### 3.1 File types

| Format | What's in it | Where to load it |
| --- | --- | --- |
| `.tdms`, raw | Group `BJRawData` with the bias (`M_E#` or `M_E1-E2`), low-gain current `M_i#a`, high-gain current `M_i#b` and piezo `U_Piezo`, plus the `Setup`, `Calibration_M#` and `Events` groups | Convert tab or Data & Traces |
| `.tdms`, with G | A `G` channel in siemens and `U_Piezo` (the `*_M.tdms` files) | Data & Traces |
| `.tdms`, other | Any channel names: type them into the channel settings (3.3) | Convert tab or Data & Traces |
| `.cvr` | Two tab-separated columns, back-calculated current (A) and piezo (nm). Old files with three columns work too | Data & Traces |
| `.txt` `.csv` `.dat` `.tsv` | Any number columns; you say which column is what | Data & Traces |
| `.pxp` (Igor packed experiment), `.ibw` (Igor binary wave) | Every numeric wave is a column, named after the wave (`folder:wave` inside data folders, `wave[0]`, `wave[1]` ... for 2D waves) | Data & Traces |

TDMS files carry the sampling rate, channel gains, piezo scale and series resistor. Text files don't, so for those you type them into the column settings (4.2). Igor waves carry their x scaling: if a wave has x in seconds, the sampling rate comes from it.

The program never recalculates what a file already gives. If there is a conductance column or channel, G is taken from it. If there is a voltage column or bias channel, it is used as the bias and the constant bias is ignored. Piezo, current and time columns are also used as they are. Only what's missing is calculated: G from the current and the bias, or the current from G and the bias.

The current can come from a single-stage amplifier (one current, Ia) or a dual-stage amplifier (a low-gain stage Ia and a high-gain stage Ib recorded at the same time). The *current amplifier* setting says which one (3.3 for TDMS, 4.2 for the other files). With a single stage there is nothing to merge and all the Ib settings are greyed out.

### 3.2 Using the tab

1. Add files or Add folder. When you click a file, the panel shows its sampling rate, length, modules, channel gains, series resistor and whether diode calibration data is in it.
2. Change the conversion settings if needed (3.3). Empty fields mean "take it from the file".
3. Preview selected shows the raw current channels, log(G/G0) with the event marks (green = approach, red = retract, grey = hold) and the piezo.
4. Then either Convert all -> .cvr, which writes one `.cvr` per TDMS file (tick Also write .csv files if you also want a table with time, G, log G, piezo, bias and current), or Load all into analysis, which converts in memory and jumps to Data & Traces.

### 3.3 Conversion settings

| Setting | Meaning | Default |
| --- | --- | --- |
| module | Which measurement module to use (`M_E#`, `M_i#a`, `M_i#b`) | 1 |
| current amplifier | `dual stage`: Ia (low gain) and Ib (high gain) are merged. `single stage`: only the Ia channel is used, nothing is merged | dual stage |
| Ia channel / Ib channel | `auto` = `M_i<module>a` / `M_i<module>b`. Type another channel name if your file uses other names; a channel that is already in A needs gain 1 | auto |
| conductance channel | `auto` uses a `G` channel only if there are no current channels. Type a channel name to take G directly from that channel, `none` to always calculate it | auto |
| conductance unit | Unit of the conductance channel: `auto` (from the channel, otherwise S), `S`, `G/G0` or `log G/G0` | auto |
| piezo channel | `auto` = `U_Piezo`; another name, or `none` | auto |
| calibration | `linear`: I = (V - offset) x gain x correction. `diode`: the piecewise polynomial stored in `Calibration_M#` (log amplifier) | linear |
| gain Ia / gain Ib (A/V) | Amplifier gains; empty = the `scale_factor` stored with the channel | from file |
| offset / correction Ia, Ib | Offset (V) and correction factor for the linear calibration | 0 / 1 |
| amplification | Extra factor on the current | 1 |
| invert current | Flips the sign | off |
| combine Ia / Ib | Dual stage only. `auto` uses the high-gain Ib unless it's saturated, then Ia. `a` or `b` uses only one | auto |
| Ib saturation (V) | Above this, Ib counts as saturated | 10 |
| Ia saturation (V) | Above this, the contact conducts more than the amplifier can measure and the sample gets the saturated G | 10.4 |
| switch to Ib below log G | Optional: only use Ib below this log(G/G0) | off |
| bias / bias channel | `measured` takes the bias channel as an array (`auto` = `M_E1-E2` if it exists, otherwise `M_E<module>`, or type a name) and the fixed bias isn't used; `fixed` uses the fixed bias. Without a bias channel the fixed bias is used | measured, auto |
| series R compensation / series R | Corrects for the current-limiting resistor, G = I / (E - I x R). Empty R = value from the Setup group | on, from file |
| saturated G (G0) | What saturated samples are set to | 100 |
| piezo scale (nm/V) | Piezo voltage to elongation; empty = the channel's `scale_factor` | from file |
| bias / current / conductance filter | none, Savitzky-Golay, median or moving average. The window is 2 x side points + 1. Can work on log values | bias: moving average, 25 points |

### 3.4 What happens in the conversion

The current is calibrated and, for a dual-stage amplifier, the two gain channels are put together (Ib where it isn't saturated, Ia elsewhere). The conductance is the current divided by the voltage that's actually across the junction, i.e. the bias minus what drops over the series resistor:

```latex
G = \frac{I}{E - I R} \cdot \frac{1}{G_0}, \qquad G_0 = \frac{2e^2}{h} = 7.748 \times 10^{-5}\ \mathrm{S}
```

If the file has a conductance channel and you choose it, none of this is done: G comes straight from the channel and only the piezo and bias are read. If the resistor takes the whole bias (E - I x R <= 0), or the low-gain amplifier saturates, the junction is in contact and beyond what can be measured, so the sample is set to the saturated G.

The `.cvr` file contains the back-calculated current, I\_bc = G x G0 x E, which is the current you'd get without the series resistor. When a `.cvr` is loaded, G is calculated back as I\_bc / bias / G0, so the bias in the column settings has to be right.

About the diode calibration: `Calibration_M#` holds five polynomials (PH, PL, M, NL, NH, coefficients in ascending order), an offset and two thresholds. With u = V - offset, the middle polynomial gives I directly when |u| is below the low threshold; outside that, the low or high polynomial gives log10 |I|, and negative u uses the negative-side polynomials with a minus sign. I haven't checked this against a measurement with module 2 yet, so test it with a known resistor before trusting it.

## 4. Data & Traces tab

Here the recordings are loaded and cut into single breaking traces (making traces too, if you want them). All the other tabs only see the traces that are accepted here.

### 4.1 Loading

Load files takes one or several files, Load folder takes every `.tdms`, `.cvr`, `.txt`, `.csv`, `.dat`, `.tsv`, `.pxp` and `.ibw` in a folder. Everything you load is analysed together. TDMS files are converted with the settings from the Convert tab (3.3), text files with the column settings below. Remove drops the selected recordings and their traces, Clear drops everything. Show columns lists the columns (or Igor waves) of a file with their numbers and lengths, so you know what to type into the column settings. The tooltip of a loaded recording says which quantities came from the file and which were calculated. The list shows how long each recording is and whether it has a piezo signal; click one to see it in the overview plot.

### 4.2 Column settings (text, .cvr and Igor files)

A column is given by its number (0 is the first) or by its name: the header of a text file or the wave name in an Igor file. Empty means the file doesn't have it.

| Setting | Meaning | Default |
| --- | --- | --- |
| conductance column | G given directly. If set, G is read from it and not calculated | empty |
| conductance unit | `G/G0`, `S` or `log G/G0` | G/G0 |
| current amplifier | `single stage`: one current column. `dual stage`: Ia and Ib, merged | single stage |
| current Ia column | The current (single stage) or the low-gain stage Ia (dual stage) | 0 |
| current Ib column | The high-gain stage Ib (dual stage only) | empty |
| Ia scale / Ib scale (A per unit) | Turns the column into A: 1 if it is in A, the amplifier gain in A/V if it is in V | 1 |
| Ib saturated from \|Ib\| | Column value from which Ib counts as saturated and Ia is used instead; empty = 99.5 % of the largest \|Ib\| in the file | empty |
| Ia saturated from \|Ia\| | Column value from which the contact is beyond the range, G = saturated G | off |
| use Ib below log G | Optional: only use Ib below this log(G/G0) | off |
| voltage column | Measured bias in V. Replaces the constant bias, and I-V analysis needs it | empty |
| constant bias (V) | Only used without a voltage column. For `.cvr` files this must be the bias of the measurement | 0.1 |
| series R (ohm) | Series resistor, G = I / (V - I x R) / G0 | 0 |
| saturated G (G0) | What saturated samples are set to | 100 |
| piezo / distance column | Piezo position or distance; empty if there is none (then distance = point number x piezo rate, 4.5) | 1 |
| piezo scale (to nm) | Factor that turns the piezo column into nm | 1 |
| time column (s) | If there is one, the sampling rate is taken from it | empty |
| sampling rate (Hz) | Used when neither a time column nor the Igor x scaling gives it | 10000 |

What is calculated depends on what you give:

| Columns given | Conductance | Current |
| --- | --- | --- |
| conductance (+ current) | from the file | from the file, or G x G0 x V |
| Ia (single stage) | I / (V - I R) / G0 | Ia |
| Ia and Ib (dual stage) | from the merged current | Ib where not saturated, Ia elsewhere |

V is the voltage column if there is one, otherwise the constant bias.

Files with one wave or column per trace (common in Igor): put a `*` in the conductance or Ia column, e.g. `G_*`. Every matching column becomes its own recording. Other columns with a `*` are paired by the same text, so with piezo column `Z_*` the wave `G_12` gets `Z_12`. Columns of different length (Igor waves) are cut to the shortest one.

Header lines are skipped, and tab, comma, semicolon or space as separator all work.

### 4.3 How the traces are found

With a piezo signal (`segmentation = piezo`) the piezo is smoothed and its speed calculated. Each point is then either moving forward, moving back, or held, which means the speed is below *hold if speed below (fraction)* times the typical speed. Very short pieces (below *min segment*) are merged with their neighbours. Every moving piece is a candidate trace and every held piece is a hold.

A moving piece counts as breaking if the conductance drops by more than a decade from start to end, and as making if it rises. If you'd rather decide by piezo direction, set *breaking direction* to `piezo+` or `piezo-`.

A lot of measurements stop the piezo in the middle of a pull to record G(t) for the noise analysis. With *join pull-hold-pull* switched on (default) the pulls before and after the stop are joined into one trace, but the hold points themselves are left out, otherwise the time spent standing still would show up as a fake peak in the histogram. The hold is kept separately for the noise tab (9).

Without a piezo signal (`segmentation = threshold`) the contacts are found where log G goes above the high limit, and each gap between two contacts is split at its lowest point into a breaking and a making trace.

### 4.4 Cutting limits and checks

All limits are log10(G/G0). Making traces are turned around and cut the same way, so every trace goes from contact to tunnelling.

| Setting | Meaning | Default |
| --- | --- | --- |
| high limit to check | The trace has to get above this (metal contact), otherwise it's rejected | 0 |
| low limit to check | The trace has to get below this, otherwise it's rejected | -4.5 |
| start at | The cut trace starts at the last point above this; empty = the high limit | empty |
| low limit to cut | The trace ends at the first point below this after the contact ... | -5.0 |
| additional length | ... plus this many points | 100 |
| reject noisy start | Throws out traces whose first *noise points* scatter more than *noise limit* (std of log G) | off, 25, 0.2 |
| max floor noise | Optional, throws out traces whose last points are noisier than this | off |
| use traces | `break`, `make` or `both` | break |

I usually put the low limits about half a decade above the noise floor. You can see the floor in the overview plot as the flat band at the bottom.

### 4.5 Distance axis

| Setting | Meaning | Default |
| --- | --- | --- |
| distance from | `piezo` uses the smoothed piezo signal, `index` uses point number x piezo rate | piezo |
| displacement ratio | Junction displacement per nm of piezo (the MCBJ attenuation factor; 1 for STM) | 1 |
| piezo rate | nm per point, for `index` and for files without piezo | 0.001 |
| zero set | Distance zero is where the trace first drops below this log G (-0.3 is 0.5 G0) | -0.3 |
| piezo smoothing | Smoothing window for direction detection and the distance axis | 10 ms |

The piezo has to be smoothed because in my files its noise (about 0.003 nm per sample) is as big as the step between samples. For an MCBJ the displacement ratio can be calibrated with the tunnelling fit in the Plateau length tab (6.3).

### 4.6 Overview and trace browser

The overview plot shows log G against time for the selected recording, with the piezo on the right axis. Breaking traces are shaded red, making traces blue and holds grey; light shading means rejected.

Below it, the trace browser shows one trace against distance with its own histogram next to it. The line above the plot gives the trace number, type, file, length and whether it's accepted (and if not, why).

| Key / button | What it does |
| --- | --- |
| Left / Right arrow, < > | Previous / next trace |
| G, Good | Mark as good and go to the next one |
| B, Bad | Mark as bad and go to the next one |
| Auto | Go back to the automatic decision |
| Accepted only | Only browse accepted traces |
| Overlay | Plot N random accepted traces on top of each other |

Manual marks are kept until you detect the traces again.

## 5. Histograms tab

The 1D histogram of log(G/G0) over all accepted traces is where the molecular conductance shows up as a peak. The 2D histogram (log G against distance) shows how far junctions get stretched at each conductance. Both are redrawn automatically whenever the traces change.

### 5.1 Settings

| Setting | Meaning | Default |
| --- | --- | --- |
| traces | All accepted traces, or one cluster from the Clustering tab | all |
| log G min / max | Conductance range of both histograms | -6.5 / 1 |
| bins / decade | Conductance resolution | 100 |
| distance min / max, distance bins | Range and resolution of the 2D histogram | -0.5 / 2 nm, 200 |
| normalization | `per trace` divides by the number of traces, so different data sets can be compared. `counts` is raw counts, `density` has area 1 | per trace |
| 2D colormap, colour max | Colour scale of the 2D histogram. The top of the scale is a percentile of the filled bins, so a few very full bins don't wash everything out | viridis, 99.5 |
| overlay clusters in 1D | If there are clusters and *all* is selected, each cluster's histogram is drawn too, scaled by its share of the traces | on |

### 5.2 Fitting the peak

Set *fit from* and *fit to* around the peak (keep away from the G0 peak and the noise floor), choose the number of Gaussians and whether you want a constant baseline, and click Fit peaks. The fit is a sum of Gaussians in log G plus the optional constant:

```latex
N(x) = c + \sum_k A_k \exp\left(-\frac{(x - \mu_k)^2}{2\sigma_k^2}\right), \qquad x = \log_{10}(G/G_0)
```

The starting values are the highest maxima in the range. For every peak you get the centre with its error, the conductance G = 10^mu in G0, the width sigma and the FWHM (2.355 sigma) in decades. The fit is drawn in red with the conductance written next to each peak.

The program remembers the last fit; the Flicker noise tab uses it for the peak window (9.3).

If a peak has a shoulder or a long tail, two Gaussians usually describe it better than forcing one. And when you compare samples, compare the per-trace histograms, not raw counts.

## 6. Plateau length tab

The plateau length is how far a junction can be pulled while the molecule is still bridging the gap. Compared with the length of the molecule it tells you whether junctions break when the molecule is fully stretched. The tab measures it in up to three conductance windows and also fits the tunnelling decay, which I use to calibrate the distance.

### 6.1 Settings

| Setting | Meaning | Default |
| --- | --- | --- |
| traces | All accepted traces or one cluster | all |
| window 1 from / to | Conductance window (log G), normally the peak +/- 1-2 sigma | -4.5 / -2.5 |
| window 2, window 3 | More windows, e.g. for a second molecular state; empty = off | off |
| length method | `span` or `count`, see below | span |
| histogram max, bins | Range and resolution of the length histogram | 2 nm, 60 |
| ignore lengths below | Leaves traces without a plateau out of the statistics | 0 nm |
| tunnelling fit from / to | log G range for the tunnelling fit (6.3) | -5 / -1 |
| expected beta (1/nm) | Decay constant you expect, for the calibration in 6.3 | 20 |

With `span`, the length runs from the first point below the upper limit to the last point above the lower limit, so short jumps out of the window don't make the plateau shorter. That's the usual stretching length. `count` multiplies the number of points inside the window by the distance step, so jumps out of the window are not counted.

### 6.2 Results

There's one length histogram per window, and the table gives the number of traces, mean, median, standard deviation, skewness, kurtosis, geometric mean, min and max for each.

Keep in mind that the real junction is longer than the measured length by the snap-back after the gold contact breaks, roughly 0.5 nm for gold. Add that before comparing with the molecule.

### 6.3 Tunnelling decay and distance calibration

For every trace log G is fitted linearly against distance in the tunnelling range, which gives beta in G \~ exp(-beta z). The right plot is the beta histogram, the median is written under the table.

For gold in vacuum beta is about 20 per nm (one decade of conductance per \~0.1 nm), in solvents it's lower; put your value in *expected beta*. The label then tells you what to multiply the *displacement ratio* (4.5) by, measured beta / expected beta. Change the ratio and detect the traces again.

## 7. Correlation tab

The correlation map shows which conductance values tend to turn up in the same trace. That way you can see things the 1D histogram averages away, for example a molecule going through two binding geometries one after the other.

For each trace r the program makes its own histogram N\_r(G) on a common set of bins, and then calculates the Pearson correlation between every pair of bins over all traces (Makk et al. 2012):

```latex
C(G_i, G_j) = \frac{\langle \delta N_r(G_i)\, \delta N_r(G_j) \rangle_r}{\sqrt{\langle \delta N_r(G_i)^2 \rangle_r \langle \delta N_r(G_j)^2 \rangle_r}}, \qquad \delta N_r = N_r - \langle N_r \rangle_r
```

| Setting | Meaning | Default |
| --- | --- | --- |
| traces | All accepted traces or one cluster | all |
| log G min / max | Range of the map | -6 / 0.5 |
| bins / decade | Coarser than in the histograms, so each trace puts enough counts in each bin | 20 |
| colour limit | Colour scale goes from minus to plus this | 0.3 |
| colormap | Diverging scale, red positive and blue negative | RdBu\_r |

Click Calculate. The diagonal is always 1, and the second plot is the mean per-trace histogram so you know where the peaks are.

How to read it: a positive region off the diagonal at (G1, G2) means traces that stay long at G1 also stay long at G2, like two plateaus that come together. Negative means the two exclude each other, e.g. two different binding geometries. Around zero there's no relation, and empty bins have no variance and stay blank.

You need a few hundred traces. With fewer the map is mostly statistical noise (of the order of 1/sqrt(number of traces)).

## 8. Clustering tab

Clustering sorts the traces into groups of similar shape without needing a reference trace, for example molecular traces vs. pure tunnelling, or two binding geometries. Once the traces are clustered, the other tabs can work on one cluster at a time.

### 8.1 What it does

Each accepted trace is first turned into a vector of fixed length. By default that's a small 2D histogram (distance x log G), flattened and normalised to sum 1, so traces of different length and sampling can be compared (same idea as Lemmer et al. 2016 and Cabosart et al. 2019). These vectors are centred and reduced with PCA, the clustering runs on the first components, and finally the clusters are numbered either by median conductance (highest first) or by size.

### 8.2 Settings

| Setting | Meaning | Default |
| --- | --- | --- |
| feature | `2D histogram` (what I'd use), `1D histogram` (log G only) or `resampled trace` (log G on a fixed distance grid) | 2D histogram |
| log G min / max, log G bins | Conductance grid of the feature | -6 / 0.5, 28 |
| distance min / max, distance bins | Distance grid of the feature | -0.2 / 1.5 nm, 28 |
| resample points | Grid size for `resampled trace` | 128 |
| PCA components | How many components the clustering uses; 0 = raw features | 10 |
| algorithm | see 8.3 | k-means |
| number of clusters | k | 3 |
| random seed | Keeps k-means, Gaussian mixture and spectral results the same from run to run | 0 |
| order clusters by | `conductance` or `size` | conductance |
| k range for scan | Which k the silhouette scan tries | 2, 8 |

### 8.3 Which algorithm

| Algorithm | Good for | Note |
| --- | --- | --- |
| k-means | A first look, compact groups of similar size | Fast, assumes round clusters |
| gaussian mixture | Groups with different spread or orientation | Elliptical clusters, slower |
| agglomerative | Hierarchical structure | Ward linkage |
| spectral | Long or curved groups | Gets slow above \~5000 traces |

### 8.4 Running it

Run clustering uses the k you set; Run + scan k also computes the silhouette score for each k in the scan range. The plots are the scree plot, the PC1-PC2 scatter coloured by cluster, the 1D histogram of each cluster, the silhouette scan (only after Run + scan k) and a 2D histogram for every cluster. The text lists the size of each cluster and the silhouette score (-1 to 1, higher means better separated). Clear clusters removes the labels again.

Don't just take the k with the highest silhouette score. With noisy data it often keeps going up with k. I take the smallest k where the cluster 2D histograms clearly look different and make physical sense; typically that's one tunnelling cluster and one or more molecular ones.

To use a cluster, pick `cluster N` in the *traces* box of Histograms, Plateau length or Correlation. When you save results, the histograms of every cluster are written automatically. Detecting the traces again removes the clusters.

## 9. Flicker noise tab

The flicker noise is used to tell through-bond transport from through-space tunnelling. The noise power NP goes with the mean conductance as NP \~ G^n. If n is close to 1, the coupling is through-bond and the current goes through the molecular orbitals. If n is close to 2, it's through-space: electrons tunnel directly between the electrodes, and since G \~ exp(-beta d), fluctuations of the gap give dG \~ G and so NP \~ G^2. Anything in between is a mix.

Two methods are implemented, the bell-shaped 2D histogram from Adak et al. (2015) and the Theil-Sen method with a stationarity test from Morris et al. (2025). They run on the same data and you always get both results.

### 9.1 What data you need

The noise is measured on G(t) traces taken while the junction is held still. Typically you pull until you're on the molecular plateau, stop the piezo for 0.1-2 s, then carry on. The Data & Traces tab finds these holds by itself (4.3), so detect the traces first.

If your data has no holds, set *segments* to `plateau windows`; then windows inside the breaking traces are used. These aren't as stationary as real holds, so be careful with the result.

Morris et al. found that n stops depending on the sampling rate above about 5 kHz for their 90 ms traces. Shorter holds probably need faster sampling.

### 9.2 Noise power

Each analysis window is either a whole hold (after trimming) or a fixed-length piece of one. The first and last few milliseconds are cut off (*trim start*, *trim end*) because the junction is still relaxing right after the piezo stops. Then the PSD of G (in G0) is calculated, with Welch's method or as a periodogram. The mean of the trace is subtracted first (*detrend* `constant`); if you don't, a trace with a step in it gives a fake 1/f^2 spectrum. Finally the PSD is integrated over the band from *integrate from* to *integrate to*:

```latex
NP = \int_{f_1}^{f_2} S_G(f)\, df \qquad [G_0^2], \qquad f_1 = 100\ \mathrm{Hz},\ f_2 = 1000\ \mathrm{Hz}
```

100-1000 Hz is the usual choice: below that there's mechanical noise, above it amplifier noise. Whatever you choose, keep it the same for all molecules you want to compare.

### 9.3 Which windows are kept

The checks below are done in this order, and every rejected window is counted under its reason in the result text.

1. Conductance range: the mean log G of the window has to be between *accept log G from* and *to*. This removes the noise floor and gold contacts.
2. Peak window: if *peak centre* and *peak sigma* are set, the mean conductance of the first and the last *end fraction* (5 %) of the window both have to be within centre +/- *n sigma* x sigma (default 2). That way only junctions that stay on the molecular plateau are used. Peak from histogram fit copies centre and sigma from the last fit in the Histograms tab. If you publish, say which n sigma you used.
3. Kurtosis (optional): throws out windows with too high excess kurtosis, i.e. with switching or spikes.
4. Drift (optional): throws out windows where the averaged log G at the start and at the end differ by more than the limit.
5. Stationarity (optional, on in the Morris preset): Augmented Dickey-Fuller test on log10|G| with a constant term, maximum lag ceil(12 (N/100)^(1/4)) and the lag chosen by AIC. The window is kept if p <= *ADF significance* (0.05), i.e. if the unit-root hypothesis is rejected. I checked my implementation against statsmodels and it gives the same statistic, lag and p-value.

With *subtract instrument floor* on, the mean PSD of the windows below *floor = windows below log G* is integrated over the same band and subtracted from every NP.

### 9.4 How n is determined

All three estimators are always calculated. *estimator for n* only decides which one is used for the bell plot and the interpretation text.

| Estimator | How | Comment |
| --- | --- | --- |
| Theil-Sen (default) | Median of the slopes between all pairs of points in log NP vs log G, which is the same as the n where Kendall's tau of (log G, log NP/G^n) is zero. Error from Sen's 95 % confidence interval | Up to 29 % outliers don't change it. This is what Morris et al. recommend |
| OLS | Least-squares slope of log NP vs log G, the same as the n where the Pearson correlation is zero | Sensitive to outliers and long tails, tends to give n too high |
| 2D Gaussian | n where the fitted 2D Gaussian (the bell) of log(NP/G^n) vs log G has zero correlation | The graphical method of Adak et al. 2015 |

The OLS error is also checked with a bootstrap (*bootstrap samples*). If there are more than 6000 windows, Theil-Sen uses a random 6000 of them, otherwise it gets slow.

### 9.5 Presets

| Preset | What it sets |
| --- | --- |
| Morris et al. 2025 (ADF + Theil-Sen) | 5 ms trimmed at both ends, one window per hold, mean-subtracted one-sided periodogram (rectangular window), 100-1000 Hz, ADF with alpha 0.05, ends within +/- 2 sigma of the peak, no drift or kurtosis cut, Theil-Sen |
| Adak et al. 2015 (2D Gaussian bell) | 100 ms windows, Welch PSD, 100-1000 Hz, 2D Gaussian estimator |

The presets don't touch your conductance range. After the Morris preset, press Peak from histogram fit so the peak window is set.

### 9.6 Plots

| Where | What |
| --- | --- |
| Top left | The segment chosen in *show segment*, accepted windows green, rejected red |
| Top middle | Mean PSD of the accepted windows, the floor PSD, the integration band and a 1/f^alpha fit |
| Top right | log NP vs log G with the OLS and Theil-Sen lines |
| Bottom left | The bell: 2D histogram of log(NP/G^n) vs log G at the chosen n, with the fitted Gaussian as contours |
| Bottom middle | Pearson r, Kendall tau and the 2D Gaussian correlation against n; each estimator is where its curve crosses zero |
| Bottom right | Residuals of the OLS and Theil-Sen fits; long tails mean outliers |

### 9.7 What to report

n with its error, which estimator, how many traces were accepted, the frequency band, the trimming, the n sigma window and whether the ADF test was on. All of this is in `*_Noise_summary.txt` (13). If OLS and Theil-Sen differ by more than their errors there are outliers in the data; I'd trust Theil-Sen and have a look at the residual plot.

As a check I ran it on simulated data with n = 1.5 and got n\_TSE = 1.52 +/- 0.02. After adding 5 % outlier windows, Theil-Sen gave 1.40 and OLS 1.33.

## 10. I-V tab

This tab is for bias sweeps. It cuts a recording into forward and backward I-V curves and makes histograms of them. For that the bias has to be recorded over time, so you need a TDMS file with a bias channel or a text or Igor file with a *voltage column* (4.2). With constant bias you just get "No bias sweeps were found".

The bias is smoothed and its turning points found (peaks and valleys at least *min sweep amplitude* apart). Every monotonic piece between two turning points is one curve, forward if the bias goes up and backward if it goes down. The current of each curve is smoothed with a Savitzky-Golay filter and dI/dV is the ratio of the filtered derivatives of I and V.

| Setting | Meaning | Default |
| --- | --- | --- |
| file | Which recording | first |
| min sweep amplitude (V) | Smallest bias swing that counts as a sweep | 0.05 |
| smoothing points, polynomial order | Savitzky-Golay window and order for the current and dI/dV | 11, 2 |
| curve must reach V below / above | Only keeps curves that cover this range; empty = off | off |
| V bins, I bins | Resolution of the 2D histograms | 100, 100 |
| log I histogram | Histogram of log10 | I |

Top left you see up to 200 curves (forward blue, backward orange) with the averaged forward and backward curves; the other three plots are 2D histograms against V of log I, of log(G/G0) with G = I/V, and of log(dI/dV / G0). The curves are saved as V/I column pairs in `IV_DATA/IV_Curves_F.txt` and `IV_Curves_B.txt` (13).

## 11. Piezo modulation tab

If a small sine is added to the piezo, the response of the conductance at that frequency gives the local decay constant:

```latex
\beta = -\frac{d \ln G}{d z}
```

For tunnelling beta is about 2 kappa, so around 20 per nm for gold in vacuum and less in solvents. On a molecular plateau it's small, because stretching the junction hardly changes G.

The modulation frequency is taken as the strongest narrow line in the piezo spectrum, unless you type it in. A line only counts if it's at least 20 times higher than the spectrum around it, so the normal piezo ramps aren't mistaken for a modulation. The recording is then split into windows of a few periods; in each window ln G and z are detrended and demodulated with a digital lock-in, and beta is minus the in-phase ratio of the two amplitudes.

| Setting | Meaning | Default |
| --- | --- | --- |
| file | Which recording | first |
| modulation frequency (Hz) | 0 = find it automatically | 0 |
| periods per window | How long the lock-in averages | 4 |
| auto search above (Hz) | Lowest frequency the detection looks at | 10 |
| log G from / to | Only windows with mean log G in this range are shown | -6 / 1 |
| displacement ratio | Junction nm per piezo nm (like in 4.5) | 1 |

The top plots are the spectra of the piezo and of log G with the modulation frequency marked, the bottom ones beta against log G for every window and the beta histogram. Without modulation in the file you get "No piezo modulation was found".

## 12. Simulation tab

The simulation makes a fake recording where the molecular conductance, plateau length and noise exponent are known. I use it to get to know the settings, to pick cutting limits and noise settings, and to check an analysis does what it should before I trust it with real data.

Every simulated cycle goes through the same trace detection as real data and looks like this: approach until the contact snaps in at 1-4 G0; pull through the atomic steps down to 1 G0, break with snap-back, and then (with *molecule probability*) a molecular plateau, otherwise straight tunnelling decay; hold, where G fluctuates with 1/f noise whose power goes as G^n; retract into the noise floor; short hold at the end.

| Setting | Meaning | Default |
| --- | --- | --- |
| number of traces | Number of cycles | 200 |
| sampling rate, speed | Hz and junction speed in nm/s | 10000, 2 |
| displacement ratio | Junction nm per piezo nm for the piezo signal | 1 |
| tunnelling decay beta | 1/nm | 10 |
| contact size range, atomic plateau length, snap-back | Contact model | 1-4 G0, 0.08 nm, 0.3 nm |
| molecule probability | Fraction of traces with a molecular plateau | 0.7 |
| molecule log G, spread | Centre and spread (decades) of the molecular conductance | -3.5, 0.25 |
| plateau length, SD, slope | Mean length and spread in nm, slope in decades/nm | 0.4, 0.15, 0.3 |
| 2nd molecule log G, fraction | A second class, to test the clustering | off, 0 |
| hold time, pull before hold, retract after hold | s, nm, nm | 1, 1, 1 |
| noise exponent n, flicker amplitude | The real n and how big the noise is | 1.5, 0.05 |
| noise floor log G, floor noise, white noise | Measurement floor and white noise | -5.5, 0.3, 0.02 |
| bias, random seed | Bias in V; the same seed gives the same data | 0.1, 1 |

Click Generate. The recording is added (or replaces what's loaded if Replace loaded data is ticked), the traces are detected with the current settings and the status bar shows the true mean molecular conductance and n. Then check: the histogram peak should be at *molecule log G*, the mean plateau length close to *plateau length* (a little off, depending on the window), the noise exponent equal to *noise exponent* within its error, and with a *2nd molecule* set the clustering should separate the two.

## 13. Saving results and settings

File -> Save results (Ctrl+S) writes everything that has been calculated into one folder, plus the figures and a complete list of the settings that were used. You choose where the folder goes and what it's called; by default it's the name of the first recording plus `_Analysis`. If the folder is already there you're asked whether to overwrite it; No makes a new folder with `_2`, `_3` and so on at the end. Every file name starts with the folder name, written `<name>` below.

The text files are tab-separated with one header line, so Origin, Excel and Python read them directly. The figures go into a `figures` subfolder (13.1).

| File | What's in it | From |
| --- | --- | --- |
| `<name>_config.json` | All settings of all tabs, the loaded files and the trace summary | always |
| `<name>_Good_Curves.txt` | All accepted traces: trace, file, z\_nm, logG\_G0 | Data & Traces |
| `<name>_all_traces.npz` | All traces including the rejected ones (NumPy file) | Data & Traces |
| `<name>_logHist.txt` | 1D histogram, logG\_G0 and counts | Histograms |
| `<name>_3D_Histogram.txt` | 2D histogram as a matrix, rows = log G, columns = distance | Histograms |
| `<name>_3D_hist_scales.txt` | Bin centres, distance\_nm and logG\_G0 | Histograms |
| `<name>_clusterN_logHist.txt`, `..._3D_Histogram.txt`, `..._3D_hist_scales.txt` | The same for each cluster | Histograms, if there are clusters |
| `<name>_peak_fit.txt` | Centre, error, G, sigma, FWHM and amplitude of each peak | Histograms, after a fit |
| `<name>_plateau length.txt` | Plateau length of every trace, one column per window | Plateau length |
| `<name>_Plateau_length_parameters.txt` | Statistics of each window and the tunnelling beta | Plateau length |
| `<name>_Correlation.txt`, `_Correlation_scales.txt` | Correlation matrix; bin centres and mean histogram | Correlation |
| `<name>_Clustering_Results.txt` | Traces per cluster, silhouette, settings | Clustering |
| `<name>_cluster_labels.txt` | Trace, cluster, PC1-PC3 | Clustering |
| `<name>_cluster_logHist.txt` | 1D histogram of each cluster | Clustering |
| `<name>_Noise_windows.txt` | Every window: segment, start, samples, G, noise power, kurtosis, drift, ADF p, accepted | Flicker noise |
| `<name>_Noise_PSD.txt` | Frequency, mean PSD, floor PSD | Flicker noise |
| `<name>_Noise_scaling.txt` | n, Pearson r, Kendall tau, 2D Gaussian rho | Flicker noise |
| `<name>_Noise_bell_hist.txt` | Bell histogram as a matrix | Flicker noise |
| `<name>_Noise_summary.txt` | n with errors for all estimators, bell parameters, counts, settings | Flicker noise |
| `<name>_IV_DATA/IV_Curves_F.txt`, `IV_Curves_B.txt` | V/I column pairs for every forward / backward curve | I-V |
| `<name>_piezo_modulation.txt` | t, logG, amplitudes, beta and phase per window | Piezo modulation |
| `figures/<name>_<plot>.pdf` and `.png` | Every plot as its own figure (13.1) | every tab |

Plateau length and Correlation results get the cluster name added (e.g. `<name>_cluster2_Correlation.txt`) if they were calculated for one cluster.

### 13.1 Figures

Every plot is saved as its own figure, sized for a journal. Nothing is combined into one image, except that 2D histograms keep their colour bar, because you can't read them without it. Multi-panel figures can be put together afterwards in Illustrator, Inkscape or similar.

|  | Default | Nature asks for |
| --- | --- | --- |
| Width | 89 mm (one column) | 89 mm single, 120 mm 1.5 columns, 183 mm double |
| Font | Arial 7 pt, tick labels and legends 6 pt | sans-serif, 5-7 pt |
| Lines | data 0.75 pt, axes and ticks 0.5 pt | at least 0.5 pt |
| Colours | Okabe-Ito palette (colour-blind safe), viridis for maps | colour-blind friendly |
| Titles | none, the file name says what it is | no panel titles |
| Files | PDF with editable text, plus a 600 dpi PNG | vector, or at least 300 dpi |

All of this can be changed under File -> Figure export settings (width, font, font size, line widths, tick direction, formats pdf/png/tiff/svg/eps, resolution). The settings are remembered and also written to `*_config.json`.

| Figure files (`figures/<name>_...`) | What they show | Tab |
| --- | --- | --- |
| `example_traces` | Random accepted traces | Data & Traces |
| `logG_histogram`, `2D_histogram`, `clusterN_logG_histogram`, `clusterN_2D_histogram` | 1D and 2D histograms of all traces and of each cluster | Histograms |
| `plateau_length_histogram`, `tunnelling_decay_histogram` | Plateau lengths and tunnelling beta | Plateau length |
| `correlation_map`, `correlation_mean_histogram` | Correlation map and mean histogram | Correlation |
| `clustering_PCA_scree`, `clustering_PCA_scatter`, `clustering_logG_histograms`, `clustering_silhouette_scan`, `clustering_clusterN_2D_histogram` | PCA and cluster results | Clustering |
| `noise_example_segment`, `noise_PSD`, `noise_scaling_fits`, `noise_bell_histogram`, `noise_correlation_vs_n`, `noise_fit_residuals` | Flicker noise | Flicker noise |
| `IV_curves`, `IV_current_histogram`, `IV_conductance_histogram`, `IV_dIdV_histogram` | I-V curves and histograms | I-V |
| `modulation_piezo_spectrum`, `modulation_logG_spectrum`, `modulation_beta_vs_G`, `modulation_beta_histogram` | Piezo modulation | Piezo modulation |

The plots on screen are drawn by the same code, so the exported figures look the same, just at journal size. If you want the whole screen view as one picture, use the save button in the plot toolbar.

### 13.2 Settings and converted files

File -> Save settings writes all settings to a JSON file and File -> Load settings reads them back; a `*_config.json` from a results folder works too, which is handy for repeating an analysis. The `.cvr` (and `.csv`) files from the Convert tab are written separately, see 3.2.

## 14. Troubleshooting

| Problem | What to do |
| --- | --- |
| All traces rejected with "never reaches high limit" | The contact doesn't get above *high limit to check*. Lower it to just under the contact level you see in the overview. If the amplifier saturates below 1 G0, check *Ia saturation* so saturated points count as contact (3.3). |
| All traces rejected with "never reaches low limit" | The noise floor is above *low limit to check*. Put both low limits about half a decade above the floor. With holds in the measurement, leave *join pull-hold-pull* on. |
| Conductance off by a constant factor | For `.cvr` files the *constant bias* in the column settings has to match the measurement; for other files check the *Ia / Ib scale* and *conductance unit*. For TDMS, look at the bias channel and series resistor in the Convert preview. |
| Distance axis is noisy or zig-zags | More *piezo smoothing*, or *distance from* `index` with the right piezo rate. |
| Distances far too long or short (MCBJ) | Calibrate the *displacement ratio* with the tunnelling fit (6.3). |
| "No hold segments were found" | No piezo stops in the recording, or no piezo signal at all. Use *segments* = `plateau windows`, or measure with holds. |
| "Only N accepted windows" | Widen *accept log G*, make *peak sigma* or *n sigma* bigger, or switch drift, kurtosis and ADF off one at a time to see which one throws the windows out. |
| "The integration band contains fewer than 2 frequency points" | The windows are too short for the band. Make *window* longer or raise *integrate from*. |
| OLS and Theil-Sen give different n | There are outliers. Use Theil-Sen, look at the residual plot, maybe switch on the ADF test. |
| "No piezo modulation was found" | There's no modulation in the file or it's too weak. Type in the *modulation frequency* if you know it. |
| "No bias sweeps were found" | The bias is constant; I-V needs sweeps. |
| Clustering changes from run to run | Keep the *random seed* fixed (it's 0 by default). |
| Window empty or tiny on a 4K screen | Set the Windows display scaling to 150-200 %. |
| "No column ... available: ..." | The column name or number isn't in the file. The message lists what is there; Show columns (4.1) does too. Numbers start at 0. |
| Igor file: "No numeric waves found" | The file has only text waves, or igor2 is missing (`py -m pip install igor2`). |
| Error message about a file | The file is broken or isn't break-junction data. Open it in a text editor or a TDMS viewer and have a look. |

You can load several files at once; they're analysed together, the trace numbers just continue and every trace remembers which file it came from. The program never changes your data files, it only reads them and writes results to the folder you pick.

To repeat an old analysis, load the `*_config.json` from its results folder with File -> Load settings, load the same data and run the same tabs.

To check that an installation works, open a command prompt in the program folder and run `GZero.exe --selftest report.txt` (or `py GZero.py --selftest report.txt`). It runs the simulation, trace detection, histogram fit, clustering, noise analysis, I-V detection and figure export without opening a window and writes the results to `report.txt`. If everything is fine the last line is `SELFTEST OK`.

## 15. References

I checked entries 1-24 against Crossref on 8 October 2026. The table at the end lists what to cite for which analysis.

### 15.1 Break-junction technique and background

1. M. A. Reed, C. Zhou, C. J. Muller, T. P. Burgin, J. M. Tour, "Conductance of a molecular junction", *Science* 278, 252-254 (1997). [doi:10.1126/science.278.5336.252](https://doi.org/10.1126/science.278.5336.252)
2. B. Xu, N. J. Tao, "Measurement of single-molecule resistance by repeated formation of molecular junctions", *Science* 301, 1221-1223 (2003). [doi:10.1126/science.1087481](https://doi.org/10.1126/science.1087481)
3. L. Venkataraman, J. E. Klare, C. Nuckolls, M. S. Hybertsen, M. L. Steigerwald, "Dependence of single-molecule junction conductance on molecular conformation", *Nature* 442, 904-907 (2006). [doi:10.1038/nature05037](https://doi.org/10.1038/nature05037)
4. N. Agrait, A. Levy Yeyati, J. M. van Ruitenbeek, "Quantum properties of atomic-sized conductors", *Physics Reports* 377, 81-279 (2003). [doi:10.1016/S0370-1573(02)00633-6](<https://doi.org/10.1016/S0370-1573(02)00633-6>)
5. T. A. Su, M. Neupane, M. L. Steigerwald, L. Venkataraman, C. Nuckolls, "Chemical principles of single-molecule electronics", *Nature Reviews Materials* 1, 16002 (2016). [doi:10.1038/natrevmats.2016.2](https://doi.org/10.1038/natrevmats.2016.2)
6. D. Xiang, X. Wang, C. Jia, T. Lee, X. Guo, "Molecular-scale electronics: from concept to function", *Chemical Reviews* 116, 4318-4440 (2016). [doi:10.1021/acs.chemrev.5b00680](https://doi.org/10.1021/acs.chemrev.5b00680)

### 15.2 Flicker noise

7. O. Adak, E. Rosenthal, J. Meisner, E. F. Andrade, A. N. Pasupathy, C. Nuckolls, M. S. Hybertsen, L. Venkataraman, "Flicker noise as a probe of electronic interaction at metal-single molecule interfaces", *Nano Letters* 15, 4143-4149 (2015). [doi:10.1021/acs.nanolett.5b01270](https://doi.org/10.1021/acs.nanolett.5b01270) Bell-shaped NP/G^n analysis; n = 1 through-bond, n = 2 through-space.
8. A. Magyarkuti, O. Adak, A. Halbritter, L. Venkataraman, "Electronic and mechanical characteristics of stacked dimer molecular junctions", *Nanoscale* 10, 3362-3368 (2018). [doi:10.1039/C7NR08354H](https://doi.org/10.1039/C7NR08354H) Determines n from the zero of the correlation, equivalent to the OLS estimator in Chapter 9.4.
9. U. Rashid, W. Bro-Jorgensen, K. B. Harilal, P. A. Sreelakshmi, R. R. Mondal, V. Chittari Pisharam, K. N. Parida, K. Geetharani, J. M. Hamill, V. Kaliginedi, "Chemistry of the Au-thiol interface through the lens of single-molecule flicker noise measurements", *J. Am. Chem. Soc.* 146, 9063-9073 (2024). [doi:10.1021/jacs.3c14079](https://doi.org/10.1021/jacs.3c14079) A recent application of noise scaling to anchor-group chemistry.
10. J. M. F. Morris, J. Potter, D. Bates, C. Wu, C. M. Robertson, S. J. Higgins, R. J. Nichols, P. J. Low, A. Vezzoli, "Digging through the (statistical) dirt: a reproducible method for single-molecule flicker noise analysis", *J. Phys. Chem. C* 129, 4097-4104 (2025). [doi:10.1021/acs.jpcc.4c07780](https://doi.org/10.1021/acs.jpcc.4c07780) ([open access](https://pmc.ncbi.nlm.nih.gov/articles/PMC11874016/)) ADF filter, Theil-Sen estimator, trimming. Data and code: [doi:10.17638/datacat.liverpool.ac.uk/2854](https://doi.org/10.17638/datacat.liverpool.ac.uk/2854).

### 15.3 Trace statistics, correlation and clustering

11. P. Makk, D. Tomaszewski, J. Martinek, Z. Balogh, S. Csonka, M. Wawrzyniak, M. Frei, L. Venkataraman, A. Halbritter, "Correlation analysis of atomic and single-molecule junction conductance", *ACS Nano* 6, 3411-3423 (2012). [doi:10.1021/nn300440f](https://doi.org/10.1021/nn300440f)
12. M. Lemmer, M. S. Inkpen, K. Kornysheva, N. J. Long, T. Albrecht, "Unsupervised vector-based classification of single-molecule charge transport data", *Nature Communications* 7, 12922 (2016). [doi:10.1038/ncomms12922](https://doi.org/10.1038/ncomms12922)
13. D. Cabosart, M. El Abbassi, D. Stefani, R. Frisenda, M. Calame, H. S. J. van der Zant, M. L. Perrin, "A reference-free clustering method for the analysis of molecular break-junction measurements", *Applied Physics Letters* 114, 143102 (2019). [doi:10.1063/1.5089198](https://doi.org/10.1063/1.5089198) See also the Comment by K. P. Hoxha, M. G. Reuter, *Appl. Phys. Lett.* 117, 136101 (2020), [doi:10.1063/5.0013531](https://doi.org/10.1063/5.0013531).
14. N. D. Bamberger, J. A. Ivie, K. N. Parida, D. V. McGrath, O. L. A. Monti, "Unsupervised segmentation-based machine learning as an advanced analysis tool for single molecule break junction data", *J. Phys. Chem. C* 124, 18302-18315 (2020). [doi:10.1021/acs.jpcc.0c03612](https://doi.org/10.1021/acs.jpcc.0c03612) Correction: *J. Phys. Chem. C* 124, 24029-24031 (2020), [doi:10.1021/acs.jpcc.0c08721](https://doi.org/10.1021/acs.jpcc.0c08721).
15. A. Magyarkuti, N. Balogh, Z. Balogh, L. Venkataraman, A. Halbritter, "Unsupervised feature recognition in single-molecule break junction data", *Nanoscale* 12, 8355-8363 (2020). [doi:10.1039/D0NR00467G](https://doi.org/10.1039/D0NR00467G)
16. Z. Balogh, G. Mezei, N. Tenk, A. Magyarkuti, A. Halbritter, "Configuration-specific insight into single-molecule conductance and noise data revealed by the principal component projection method", *J. Phys. Chem. Lett.* 14, 5109-5118 (2023). [doi:10.1021/acs.jpclett.3c00677](https://doi.org/10.1021/acs.jpclett.3c00677) Combines PCA classes with noise analysis.

### 15.4 Statistical and signal-processing methods

17. H. Theil, "A rank-invariant method of linear and polynomial regression analysis" (originally *Proc. Koninklijke Nederlandse Akademie van Wetenschappen*, 1950); reprinted in *Henri Theil's Contributions to Economics and Econometrics*, Advanced Studies in Theoretical and Applied Econometrics, pp. 345-381 (Springer, 1992). [doi:10.1007/978-94-011-2546-8\_20](https://doi.org/10.1007/978-94-011-2546-8_20)
18. P. K. Sen, "Estimates of the regression coefficient based on Kendall's tau", *J. Am. Stat. Assoc.* 63, 1379-1389 (1968). [doi:10.1080/01621459.1968.10480934](https://doi.org/10.1080/01621459.1968.10480934)
19. M. G. Kendall, "A new measure of rank correlation", *Biometrika* 30, 81-93 (1938). [doi:10.1093/biomet/30.1-2.81](https://doi.org/10.1093/biomet/30.1-2.81)
20. D. A. Dickey, W. A. Fuller, "Distribution of the estimators for autoregressive time series with a unit root", *J. Am. Stat. Assoc.* 74, 427-431 (1979). [doi:10.1080/01621459.1979.10482531](https://doi.org/10.1080/01621459.1979.10482531)
21. J. G. MacKinnon, "Approximate asymptotic distribution functions for unit-root and cointegration tests", *J. Business & Economic Statistics* 12, 167-176 (1994). [doi:10.1080/07350015.1994.10510005](https://doi.org/10.1080/07350015.1994.10510005)
22. P. Welch, "The use of fast Fourier transform for the estimation of power spectra: a method based on time averaging over short, modified periodograms", *IEEE Trans. Audio Electroacoust.* 15, 70-73 (1967). [doi:10.1109/TAU.1967.1161901](https://doi.org/10.1109/TAU.1967.1161901)
23. A. Savitzky, M. J. E. Golay, "Smoothing and differentiation of data by simplified least squares procedures", *Analytical Chemistry* 36, 1627-1639 (1964). [doi:10.1021/ac60214a047](https://doi.org/10.1021/ac60214a047)

### 15.5 Software libraries

24. P. Virtanen et al., "SciPy 1.0: fundamental algorithms for scientific computing in Python", *Nature Methods* 17, 261-272 (2020). [doi:10.1038/s41592-019-0686-2](https://doi.org/10.1038/s41592-019-0686-2)
25. F. Pedregosa et al., "Scikit-learn: machine learning in Python", *J. Machine Learning Research* 12, 2825-2830 (2011). No DOI; [jmlr.org](https://jmlr.org/papers/v12/pedregosa11a.html). Not in Crossref and not re-checked.

### 15.6 What to cite for each analysis

| Analysis in this program | Cite |
| --- | --- |
| Conductance histograms, plateau length | 2, 3 (technique); 4 or 5 (review) |
| 2D cross-correlation | 11 |
| Clustering | 12, 13; 15 or 14 for related methods |
| Flicker noise, bell / 2D Gaussian | 7 |
| Flicker noise, OLS / zero correlation | 8 |
| Flicker noise, ADF + Theil-Sen | 10, with 17, 18, 20 |
| PSD (Welch), smoothing | 22, 23 |
| Libraries used | 24, 25 |
