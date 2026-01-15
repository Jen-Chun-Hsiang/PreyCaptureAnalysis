# Methods (analysis pipeline)

This document describes the *conceptual* and *mathematical* analysis steps used to transform raw electrophysiology recordings and visual stimulus descriptions into receptive-field estimates and model-based predictions of RGC responses. The intent is that an independent reader could reproduce the analysis from scratch using only the raw data and stimulus timing.

## Overview

We analyzed retinal ganglion cell (RGC) responses to several classes of visual stimuli:

- **Moving noise (jittered checkerboard white noise)** used to estimate linear spatiotemporal receptive fields and static nonlinearities.
- **Spot stimuli (multiple diameters)** used to quantify center–surround organization.
- **Moving bar stimuli** used to test whether receptive-field parameters estimated from noise generalize to structured motion, and to fit reduced models including optional surround drive and adaptive/divisive components.

Analyses were performed per recording/cell and then aggregated across groups defined by response polarity (ON vs OFF) and retinal location (e.g., nasal vs temporal).

### Moving-noise stimulus (jittered checkerboard)

The “moving noise” stimulus is an extension of a standard checkerboard white-noise stimulus. Rather than presenting an independent checkerboard on a fixed pixel grid, the entire noise canvas is **slightly jittered (sub-tile shifts) each frame**. This produces finer effective spatial sampling while still using relatively large checkerboard elements that reliably drive spikes. Conceptually, this is equivalent to dithering the stimulus grid over time so that spike-triggered averaging can recover more detailed receptive-field structure than would be possible with a static grid at the same checker size.

## Spike extraction and response representation

### Spike detection

Raw voltage traces were detrended and spikes were detected using a peak/threshold-based method. Detection thresholds were chosen per recording (to account for variability in noise floor and spike amplitude). Detected spike times were converted to a binned spike train.

### Firing-rate representations (PSTHs): binning and smoothing parameters

Across the stimulus-processing scripts (moving bars, moving spots, and stationary/varied-size spots), firing-rate traces were computed from detected spike times using **fixed binning** followed by an **optional smoothing** step:

- **Raw sampling rate**: electrophysiology traces were sampled at $10\,\mathrm{kHz}$.
- **Fixed binning**: spike counts were computed in non-overlapping bins on a grid of $F_s = 100\,\mathrm{Hz}$ (bin width $\Delta t = 10\,\mathrm{ms}$) by summing spikes within each $10\,\mathrm{ms}$ bin.
- **Gaussian smoothing for PSTHs**: the binned spike counts were smoothed using a Gaussian window with width **$50\,\mathrm{ms}$** (MATLAB `smoothdata(..., 'gaussian', 0.05, 'SamplePoints', t)`), and then converted to firing rate in spikes/s by multiplying by $F_s$.

This smoothing was used to produce visually stable PSTHs and to define a continuous-valued firing-rate signal for downstream analyses on these stimuli.

### Time base and downsampling

To facilitate stimulus–response modeling, spike trains and stimulus time series were represented on a common discrete time grid. In the current analysis, responses were commonly downsampled/binned to $F_s = 100\ \mathrm{Hz}$ (bin width $\Delta t = 10\ \mathrm{ms}$). All temporal filtering and model fitting were performed on this grid unless otherwise noted.

### Repeat averaging and mean firing rate

Many stimulus conditions were repeated multiple times. For each condition, we computed:

- A **trial-averaged PSTH** by averaging the per-trial firing-rate traces across repeats at each time bin.
- A **mean firing rate** within a defined analysis window as the mean of the trial-averaged PSTH over time (equivalently, the average over both time and repeat when the window is the same for all repeats).

Unless otherwise stated for a specific analysis, reported “mean firing rate” values refer to this repeat-averaged quantity (to reduce trial-to-trial variability).

### Alignment to stimulus timing

Stimulus frames (and/or stimulus parameter time series) were aligned to the response time base via experiment-specific synchronization (stimulus timestamps and acquisition timing). When necessary, stimulus signals were resampled to match the response binning.

## Moving-noise (white-noise) receptive-field estimation

### Spatiotemporal stimulus representation

Let $S(\mathbf{x}, t)$ denote the stimulus value at pixel/location $\mathbf{x}$ at time bin $t$ (e.g., intensity or contrast). The stimulus history window preceding a spike spans a fixed number of bins corresponding to a time window (e.g., hundreds of ms).

**Implementation details (code-based parameters).** For the moving-noise RF estimation, stimuli and responses were represented on a $F_s=100\,\mathrm{Hz}$ time grid (bin width $\Delta t=10\,\mathrm{ms}$), and the spike-triggered history window was

$$
\mathrm{WinT} = [-0.5,\ 0]~\mathrm{s},
$$

corresponding to $N_T=\lceil 0.5\,F_s\rceil = 50$ time bins.

The moving-noise stimulus frames were reconstructed from the experiment seed using MATLAB’s `RandStream('mrg32k3a','seed', IN.NoiseUniqueSeed)` to regenerate the underlying binary checkerboard, then:

- Upsampled to a fixed image grid of $800\times 600$ pixels using nearest-neighbor interpolation.
- Translated each frame by the recorded per-frame jitter offsets (`OUT.MovingSteps`) by shifting indices (out-of-bounds pixels filled with 0).
- Converted from binary values to contrast values in $\{-1, +1\}$ via $S \leftarrow 2(S-0.5)$.

### Spike-triggered average (STA)

The spatiotemporal spike-triggered average estimates the mean stimulus preceding spikes:

$$
\mathrm{STA}(\mathbf{x}, \tau) = \frac{1}{N_{\mathrm{sp}}} \sum_{k=1}^{N_{\mathrm{sp}}} S\big(\mathbf{x},\ t_k - \tau\big),
$$

where $t_k$ are spike times on the binned grid, $\tau \ge 0$ indexes time lag, and $N_{\mathrm{sp}}$ is the number of spikes used.

For spatial summaries (e.g., for center localization), the STA was reduced across time lags (e.g., by selecting or integrating over informative lags) and normalized to yield a spatial STA “strength” image.

**Implementation details (code-based parameters).** Spikes were detected on the $10\,\mathrm{kHz}$ voltage trace using a peak-finding step (MATLAB `findpeaks` on the inverted trace) with a refractory constraint `MinPeakDistance = 40` samples ($4\,\mathrm{ms}$ at $10\,\mathrm{kHz}$). Spike times were then binned to the $F_s=100\,\mathrm{Hz}$ grid by summing spikes in non-overlapping $10\,\mathrm{ms}$ bins.

The STA was accumulated using **per-bin spike counts** (not just spike times). For each $10\,\mathrm{ms}$ bin $t$ with spike count $r(t)$, the code adds $r(t)$ copies of the preceding stimulus history cube spanning the last 50 bins (covering $0.5\,\mathrm{s}$):

$$
\mathrm{STA}(\mathbf{x},\tau) \propto \sum_{t} r(t)\,S(\mathbf{x}, t-\tau), \quad \tau\in\{1,\dots,50\}.
$$

Within each stimulus block, the accumulated STA cube was normalized by the total spike count in that block (with `+eps` for numerical stability) and then averaged across noise blocks.

To avoid circularity when estimating the static nonlinearity, the code randomly withheld a fraction of response bins from the STA computation (`nonlinear_ratio = 0.25`, i.e., 25% of eligible bins), reserving those withheld bins for the nonlinearity estimate.

### Spatial RF localization and parametric Gaussian fit

To obtain compact receptive-field (RF) descriptors, the spatial RF was fit with a 2D Gaussian model (optionally elliptical and rotated). In general form:

$$
G(\mathbf{x}) = A\,\exp\Big(-\tfrac{1}{2}(\mathbf{x}-\boldsymbol\mu)^\top\Sigma^{-1}(\mathbf{x}-\boldsymbol\mu)\Big) + b,
$$

with amplitude $A$, center $\boldsymbol\mu$, covariance $\Sigma$ (controlling size and ellipticity), and baseline $b$. Parameters were estimated by minimizing a squared-error objective over pixels inside a mask (masking excludes invalid or low-confidence regions):

$$
\hat\theta = \arg\min_{\theta} \sum_{\mathbf{x}\in\mathcal{M}} \big(I(\mathbf{x}) - G_\theta(\mathbf{x})\big)^2,
$$

where $I(\mathbf{x})$ is the spatial STA-derived image and $\mathcal{M}$ is the spatial mask.

From $\Sigma$ we derived size-related quantities (e.g., equivalent diameter estimates) and ellipticity (ratio of principal-axis standard deviations).

**Implementation details (code-based parameters).** In the moving-noise scripts, the spatial RF “strength” image was computed as the per-pixel standard deviation across STA lags,

$$
I(\mathbf{x}) = \mathrm{std}_{\tau}(\mathrm{STA}(\mathbf{x},\tau)).
$$

This image was median-filtered (`medfilt2`) and then thresholded to form a spatial mask $\mathcal{M}$:

- A helper `peak_distance` was applied to the filtered image.
- Only peaks beyond a minimum distance `minD = 100` (pixels) were considered.
- A threshold was set at the 99th percentile of the corresponding peak scores: `ythr = quantile(y(x>minD), 0.99)`.
- A binary mask was formed as `binary_image = smtstdSTA > ythr`, and the **largest 4-connected component** was kept via `largest_segment_4conn_mask`.

Spatial RF parameters (center location, size, and orientation) were obtained by fitting the spatial image with a rotated elliptical Gaussian model implemented in `gaussian2d.m`. The fitted model is a **sum of two rotated elliptical Gaussians plus a constant offset**:

$$
G(x,y)= A\exp\Big(-\tfrac{1}{2}\big(\tfrac{x_\theta^2}{\sigma_x^2}+\tfrac{y_\theta^2}{\sigma_y^2}\big)\Big)
\; +\; A_s\exp\Big(-\tfrac{1}{2}\big(\tfrac{x_\theta^2}{\sigma_{x,s}^2}+\tfrac{y_\theta^2}{\sigma_{y,s}^2}\big)\Big)
\; +\; b,
$$

where $(x_\theta, y_\theta)$ are coordinates rotated by angle $\theta$. The parameter vector is

$$
	heta_{\mathrm{spatial}} = [x_0, y_0, \sigma_x, \sigma_y, \theta, b, A, \sigma_{x,s}, \sigma_{y,s}, A_s].
$$

Fits used Nelder–Mead optimization (`fminsearch`) with objective

$$
\hat\theta = \arg\min_{\theta}\; \big(1 - \mathrm{corr}(I(:), G_{\theta}(:))\big),
$$

starting from the initialization

$$
[x_0,y_0,\sigma_x,\sigma_y,\theta,b,A,\sigma_{x,s},\sigma_{y,s},A_s]
= [\tfrac{W}{2},\tfrac{H}{2},50,50,0,0.1,1,200,200,0.1]
$$

(in pixels, with a small random jitter added to $(x_0,y_0)$), and with `MaxFunEvals = 600 * nParams`.

Reported spatial RF size and orientation were derived from the fitted center Gaussian parameters $(\sigma_x,\sigma_y,\theta)$ and converted to physical units using the OLED calibration (`OLED.pixelSize`, i.e., $\mu\mathrm{m}$/pixel).

### Temporal receptive field (temporal filter)

A temporal filter $h(\tau)$ was computed by spatially pooling the spatiotemporal STA with the fitted/selected spatial RF weights. Conceptually:

$$
\mathrm{tRF}(\tau) = \sum_{\mathbf{x}\in\mathcal{M}} w(\mathbf{x})\,\mathrm{STA}(\mathbf{x},\tau),
$$

where $w(\mathbf{x})$ are spatial weights (e.g., masked spatial RF values).

To provide a smooth, low-dimensional description, temporal filters were fit using a **difference-of-Gaussians in time**, i.e., a biphasic form:

$$
\hat h(\tau) = a_1\exp\Big(-\frac{(\tau-\mu_1)^2}{2\sigma_1^2}\Big)
- a_2\exp\Big(-\frac{(\tau-\mu_2)^2}{2\sigma_2^2}\Big) + b.
$$

Fit parameters include peak times ($\mu_1,\mu_2$), widths ($\sigma_1,\sigma_2$), amplitudes ($a_1,a_2$), and offset $b$. From $\hat h(\tau)$ we computed standard temporal metrics such as time-to-peak, temporal width, and degree of biphasy.

**Implementation details (code-based parameters).** The temporal filter was computed by (i) median filtering the STA cube in space–time (`medfilt3`), then (ii) projecting the filtered STA cube onto the masked spatial RF weights (the masked, median-filtered $I(\mathbf{x})$):

$$
\mathrm{tRF}(\tau) = \sum_{\mathbf{x}} \mathrm{STA}_{\mathrm{filt}}(\mathbf{x},\tau)\,w(\mathbf{x}),
$$

with $\tau$ sampled on the same $F_s=100\,\mathrm{Hz}$ grid over $[-0.5,0)$ seconds.

Temporal fits were performed with the helper `GaussianTemporalFilter.m`, which normalizes the empirical temporal filter by its maximum (`scaling = max(tRF)`, storing that scaling factor) and then fits a two-Gaussian biphasic form by constrained optimization (`fmincon`). The fitted parameter vector is

$$
w = [\sigma_1,\sigma_2,\mu_1,\mu_2,a_1,a_2,b],
$$

with initial values $[2,\ 2,\ \arg\max(tRF),\ \arg\min(tRF),\ 0,\ 0,\ 0]$ (in **samples**) and box constraints:

- $\sigma_1,\sigma_2\in[0,20]$ samples
- $\mu_1,\mu_2\in[1,N_T]$ samples
- $a_1,a_2\in[0,10]$
- $b\in[-1,1]$

To report time-to-peak with finer resolution, the downstream summary code interpolates $\mathrm{tRF}$ to 1000 samples over the $[-0.5,0]$ s window (cubic interpolation) and takes the ON peak as the maximum (OFF peak as the minimum), reporting peak time in ms. Temporal “width” is computed as the total duration where the interpolated trace exceeds a fixed threshold magnitude (`hwith_thr = 0.5`, with sign depending on ON/OFF).

#### Biphasic index (code-based definitions)

In the current analysis code, biphasy of the temporal filter was quantified in two closely related ways.

1) **Peak-balance biphasic index (from the empirical temporal filter)**. Let $h(t)$ denote the empirical temporal filter (the spatially pooled STA-derived temporal receptive field, `tRF`). The code first interpolates $h(t)$ to a finer temporal grid (cubic interpolation to 1000 samples over the lag window) and then computes the positive and negative peak magnitudes:

$$
p = \max_t h(t), \qquad n = \left|\min_t h(t)\right|.
$$

The peak-balance biphasic index is then

$$
\mathrm{BI}_{\mathrm{peaks}} = 1 - \frac{|p-n|}{p+n} 
= \frac{2\,\min(p,n)}{p+n},
$$

which lies in $[0,1]$ (0 = purely monophasic, 1 = equal-magnitude positive and negative lobes).

2) **Fitted-amplitude biphasic index (from the two-Gaussian fit)**. For the fitted model

$$
\hat h(t) = a_1\,\mathcal{N}(t;\mu_1,\sigma_1) - a_2\,\mathcal{N}(t;\mu_2,\sigma_2) + b,
$$

the analysis defines a biphasic “strength” ratio from the fitted amplitudes. There are two closely related parameterizations in the codebase:

- **Independent-amplitude parameterization (`GaussianTemporalFilter.m` + `gaussian_temporalfilter.m`)**. The fitted parameter vector is

$$
w=[\sigma_1,\sigma_2,\mu_1,\mu_2,a_1,a_2,b],
$$

and the fitted filter is

$$
\hat h(t)= a_1\,\mathcal{N}(t;\mu_1,\sigma_1) - a_2\,\mathcal{N}(t;\mu_2,\sigma_2) + b.
$$

- **Amplitude+ratio parameterization (`GaussianTemporalFilter2.m`, used by the white-noise summary script)**. The fitted parameter vector is

$$
w=[\sigma_1,\sigma_2,\mu_1,\mu_2,A,\rho,b],
$$

and the fitted filter is

$$
\hat h(t)= A\,\mathcal{N}(t;\mu_1,\sigma_1) - (\rho A)\,\mathcal{N}(t;\mu_2,\sigma_2) + b.
$$

Here $A$ is the primary-lobe amplitude scale and $\rho$ is the secondary-to-primary lobe amplitude ratio (constrained to be nonnegative in the optimizer).

For paper reporting, the natural “ratio of secondary to primary lobe magnitude” from the fitted model is

$$
\mathrm{BI}_{\mathrm{strength}} = \frac{\min(|a_1|,|a_2|)}{\max(|a_1|,|a_2|)}\in[0,1],
$$

which reduces to $\min(\rho, 1/\rho)$ under the amplitude+ratio parameterization.

**As implemented in the current analysis scripts**, the fitted-strength biphasy metric is computed from the stored optimizer output `Gauss_TF_est` as a ratio of columns 5 and 6 with an ON/OFF-specific inversion:

- For ON cells: `TF_biphasic_stregth = Gauss_TF_est(:,6) ./ Gauss_TF_est(:,5)`
- For OFF cells: `TF_biphasic_stregth = Gauss_TF_est(:,5) ./ Gauss_TF_est(:,6)`

This choice enforces a consistent ordering across polarities (so values remain $\le 1$ if columns 5 and 6 represent the two lobe magnitudes). When using `GaussianTemporalFilter2.m`, note that column 6 is the ratio $\rho$ and the inhibitory-lobe amplitude is $\rho\,A$.

$$
\mathrm{BI}_{\mathrm{strength}} \approx \frac{\text{secondary lobe magnitude}}{\text{primary lobe magnitude}}.
$$

## Static nonlinearity (LN output stage)

### Generator signal

Using the estimated linear spatiotemporal RF, we computed a scalar **generator signal** $g(t)$ by projecting the stimulus onto the filter:

$$
 g(t) = \sum_{\tau=0}^{T_h} \sum_{\mathbf{x}} K(\mathbf{x})\,h(\tau)\,S(\mathbf{x}, t-\tau),
$$

where $K(\mathbf{x})$ is a spatial RF (e.g., Gaussian fit or masked RF) and $h(\tau)$ is the temporal filter.

### Nonlinearity estimation

The static nonlinearity maps generator signal to firing rate. We estimated it empirically by binning $g(t)$ values and computing the mean firing rate within each bin. Here $r(t)$ denotes the measured firing rate on the discrete time grid (typically binned at $100\,\mathrm{Hz}$; smoothing may be applied depending on the stimulus/analysis stage):

$$
\hat r(g_i) = \frac{1}{|\mathcal{T}_i|}\sum_{t\in\mathcal{T}_i} r(t), \quad \mathcal{T}_i = \{t : g(t)\in\text{bin }i\}.
$$

We then fit a smooth parametric curve to the binned estimate. A common choice in this pipeline is a scaled cumulative normal (sigmoidal) form:

$$
 r(g) = A\,\Phi\Big(\frac{g-\mu}{\sigma}\Big) + b,
$$

where $\Phi$ is the standard normal CDF, and $(A,\mu,\sigma,b)$ are fitted parameters. For OFF cells, sign conventions can be handled by flipping the generator signal or the fitted curve so that the fitted nonlinearity is monotonically increasing in the effective drive.

Derived quantities (used for comparisons) include threshold-like measures (e.g., an effective $\mu$) and slope/softness measures (e.g., related to $\sigma$).

## Spot-size / center–surround analysis

Spot-size stimuli probe spatial integration and surround suppression.

### ON/OFF cell classification

Cell polarity (ON vs OFF) was determined during patch-clamp experiments using a step-like spot stimulus (spot turns ON for a fixed epoch and then turns OFF), presented across multiple spot sizes. Cells were labeled:

- **ON** if their response was dominated by the **ON (light increment) phase** of the spot.
- **OFF** if their response was dominated by the **OFF (light decrement) phase** of the spot.

This classification was used for grouping and comparisons throughout the downstream analyses.

### Per-size response extraction

For each spot diameter $d$, we computed a response statistic (typically mean firing rate) during a fixed stimulus window. This yields a size-tuning curve $R(d)$ per cell.

### Size Index (SI)

A simple center–surround suppression metric was computed from the size-tuning curve using two regimes:

- “Center” sizes: smaller spots below a threshold (e.g., $d<800\ \mu\mathrm{m}$).
- “Large/surround” sizes: large spots (e.g., $d\in\{800,1200\}\ \mu\mathrm{m}$).

Let $C$ be an estimate of the center-dominated response (computed from the peak in the small-size regime and an adjacent size), and let $S$ be the mean response at the large sizes. The size index is:

$$
\mathrm{SI} = \frac{S - C}{C}.
$$

Negative SI indicates suppression for large spots relative to the center response.

### Parametric DoG-like fits of size tuning

To obtain interpretable center and surround scale parameters, size-tuning curves were fit with a difference-of-Gaussian–motivated integrated form (equivalent to integrating a 2D Gaussian over a disk/diameter). One practical parameterization used is:

$$
\hat R(d) = b
+ k_c\Big(1 - \exp\big(-\tfrac{((d-\mu)/2)^2}{2\sigma_c^2}\big)\Big)
- k_s\Big(1 - \exp\big(-\tfrac{((d-\mu)/2)^2}{2\sigma_s^2}\big)\Big),
$$

where $(k_c,\sigma_c)$ describe center gain/scale, $(k_s,\sigma_s)$ describe surround gain/scale, $b$ is baseline, and $\mu$ is an optional diameter offset. In some fits, the ratio $\sigma_s/\sigma_c$ is constrained (e.g., fixed to a constant) to stabilize estimation.

A derived surround-to-center weight ratio is computed from the fitted parameters (proportional to integrated weights):

$$
\mathrm{WR} = \frac{k_s\,\sigma_s^2}{k_c\,\sigma_c^2}.
$$

### Bootstrap estimation of RF size parameters

For robust estimation and confidence intervals, center and surround size parameters were also estimated via bootstrap resampling of the size-tuning data. Each bootstrap draw refits the model and yields distributions over $(\sigma_c,\sigma_s)$ and derived size measures (e.g., diameters at a fixed fraction of saturation). Confidence intervals are reported as percentiles of the bootstrap distribution.

## Moving-bar modeling

Moving bars provide a structured stimulus for testing generalization of RF-based models.

### Linear drive from noise-derived RFs

We constructed a spatiotemporal linear filter from the noise-derived spatial RF and temporal filter, then applied it to the moving-bar stimulus to compute a predicted linear drive.

We explicitly computed two drives:

- **Center drive** $x(t)$: projection using a “center” spatial filter.
- **Surround drive** $x_s(t)$: projection using a broader spatial filter representing surround integration.

Conceptually:

$$
 x(t) = \sum_{\tau}\sum_{\mathbf{x}} K_c(\mathbf{x})\,h(\tau)\,S(\mathbf{x},t-\tau),
 \qquad
 x_s(t) = \sum_{\tau}\sum_{\mathbf{x}} K_s(\mathbf{x})\,h(\tau)\,S(\mathbf{x},t-\tau).
$$

Here $K_c$ and $K_s$ are center/surround spatial filters (e.g., Gaussians with different scales; $K_s$ may be derived from $K_c$ by a scale factor).

### Rectified linear output model (LN with optional surround)

To map linear drive(s) to firing rate, we fit a rectified affine model. Center-only:

$$
\hat r(t) = \max\big(\alpha\,x(t) + \beta,\ 0\big).
$$

Center + surround suppression:

$$
\hat r(t) = \max\big(\alpha\,(x(t) - \gamma\,x_s(t)) + \beta,\ 0\big), \quad \gamma\ge 0.
$$

Parameters $(\alpha,\beta,\gamma)$ were fit by minimizing mean-squared error between predicted and measured binned firing rates over all stimulus conditions and time samples.

#### Optional CSR constraint

When an independent center–surround measure is available from spot-size analysis, we optionally regularize the surround strength toward a target value $\gamma_0$:

$$
\mathcal{L} = \frac{1}{T}\sum_t\big(\hat r(t) - r(t)\big)^2 + \lambda\,(\gamma-\gamma_0)^2,
$$

where $\lambda$ controls constraint strength. (In practice, sign conventions are handled so that inhibitory surround corresponds to positive $\gamma$.)

### Adaptive/divisive LNK-family models (optional)

To capture adaptation and gain control beyond static LN models, we fit reduced kinetic/divisive models with an internal state and divisive normalization.

A representative 1-state model is:

State update (adaptation):

$$
 a_{t+1} = a_t + \Delta t\,\frac{\alpha_d\,F(x_t) - a_t}{\tau},
 \qquad F(x)=\max(0, x-\theta).
$$

Divisive normalization and combined drive:

$$
 \mathrm{den}_t = \sigma_0 + \alpha\,a_t,
$$

$$
 \tilde y_t = \frac{x_t}{\mathrm{den}_t} + w_{xs}\,\frac{x_{s,t}}{\mathrm{den}_t} + \beta\,a_t + b_{out}.
$$

Output nonlinearity (default softplus):

$$
 \hat r_t = \mathrm{softplus}(g_{out}\,\tilde y_t),
 \qquad \mathrm{softplus}(z)=\log(1+e^z).
$$

Parameters $(\tau,\alpha_d,\theta,\sigma_0,\alpha,\beta,b_{out},g_{out},w_{xs})$ were optimized to minimize prediction error, optionally with a soft regularizer that encourages $w_{xs}$ to match a target derived from center–surround measurements.

### Model evaluation and repeat reliability

Because responses were measured with repeated presentations, we used **correlation-based repeat reliability** as a quality-control metric. The goal is to quantify how much of the time-varying response is stimulus-driven (repeatable) versus dominated by intrinsic trial-to-trial variability.

In all cases below, reliability is based on the Pearson correlation coefficient computed between firing-rate time series from repeated presentations of the *same* stimulus condition.

#### Moving bar and moving spot (condition-wise response quality)

For moving bar and moving spot stimuli, responses were organized into a tensor
repeats $\times$ time for each condition (direction/speed/size, etc.). For a given condition we assembled the per-repeat PSTHs into a matrix
$\mathbf{Y}\in\mathbb{R}^{n_{\mathrm{rep}}\times T}$ (after removing time bins that were NaN for any repeat due to variable trial lengths). We then computed the repeat-by-repeat correlation matrix:

$$
\mathbf{R} = \mathrm{corr}(\mathbf{Y}^\top),
$$

and defined a scalar repeatability score as the mean of the off-diagonal entries (i.e., the mean pairwise correlation across repeats):

$$
\rho = \frac{2}{n_{\mathrm{rep}}(n_{\mathrm{rep}}-1)}\sum_{i<j} R_{ij}.
$$

In the current processing scripts this is typically stored/visualized as a squared correlation ("$R^2$"):

$$
\rho_{R^2} = \frac{2}{n_{\mathrm{rep}}(n_{\mathrm{rep}}-1)}\sum_{i<j} R_{ij}^2,
$$

and plotted as $\log(\rho_{R^2})$ as a compact "response quality" summary across stimulus conditions.

#### Moving noise (repeated-block reliability)

For moving-noise stimuli, the stimulus sequence includes explicitly repeated blocks. Repeat reliability was computed by extracting the response to each repeated block, binning spikes at $F_s = 100\,\mathrm{Hz}$ (bin width $10\,\mathrm{ms}$) into a per-block response trace, optionally smoothing each trace, and then computing the mean pairwise correlation across repeated blocks:

$$
\rho_{\mathrm{WN}} = \frac{2}{n_{\mathrm{rep}}(n_{\mathrm{rep}}-1)}\sum_{i<j} \mathrm{corr}(y_i, y_j).
$$

In the current moving-noise repeatability script, smoothing is implemented by convolving the binned spike-count trace with a Gaussian kernel (`gaussian_smooth_1d`) using a kernel length of 100 samples and Gaussian width parameter $\sigma=0.08\,\mathrm{s}$ (when enabled).

#### Stationary/varied-size spot stimuli

For stationary/varied-size spots, the main analysis scripts primarily use repeats to estimate mean response and variability (e.g., SEM) of firing rates within fixed ON/OFF analysis windows. A correlation-based repeat-reliability metric can be computed analogously to the moving-bar/moving-spot procedure (treating each repeat’s PSTH as a vector and averaging pairwise correlations), but this is not the primary QC summary implemented in the current stationary-spot processing script.

#### Model prediction vs. repeatable response (repeat-aware correlation)

For moving-bar model evaluation, prediction performance is computed using a repeat-aware correlation: predictions are compared to the measured response while respecting repeat structure by computing correlations within each repeat and then averaging these correlations across repeats (omitting NaNs).

Performance was summarized across stimulus conditions and then aggregated across cells/groups.

## Grouping and aggregation

Cells were grouped by response polarity (ON vs OFF; defined from the step-spot ON vs OFF phase response during patching) and by retinal location categories (e.g., nasal vs temporal). For each fitted quantity (spatial RF size/ellipticity, temporal filter parameters, nonlinearity parameters, center–surround indices, moving-bar model parameters and performance), we computed group means and uncertainty estimates (e.g., SEM across cells), and compared groups using appropriate statistical tests when needed.

## Practical reproducibility notes

- The pipeline relies on consistent stimulus–response alignment and careful spike detection; both materially affect STA estimates and downstream model fits.
- Regularization/constraints (e.g., fixed surround-to-center scale ratios for spot-size DoG fits; CSR penalties for moving-bar fits) are used to stabilize estimation given finite data.
- Wherever bootstrap resampling is used, confidence intervals should be reported from the bootstrap distribution (percentiles).
