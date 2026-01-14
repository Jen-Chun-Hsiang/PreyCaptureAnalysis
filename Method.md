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

### Time base and downsampling

To facilitate stimulus–response modeling, spike trains and stimulus time series were represented on a common discrete time grid. In the current analysis, responses were commonly downsampled/binned to $F_s = 100\ \mathrm{Hz}$ (bin width $\Delta t = 10\ \mathrm{ms}$). All temporal filtering and model fitting were performed on this grid unless otherwise noted.

### Alignment to stimulus timing

Stimulus frames (and/or stimulus parameter time series) were aligned to the response time base via experiment-specific synchronization (stimulus timestamps and acquisition timing). When necessary, stimulus signals were resampled to match the response binning.

## Moving-noise (white-noise) receptive-field estimation

### Spatiotemporal stimulus representation

Let $S(\mathbf{x}, t)$ denote the stimulus value at pixel/location $\mathbf{x}$ at time bin $t$ (e.g., intensity or contrast). The stimulus history window preceding a spike spans a fixed number of bins corresponding to a time window (e.g., hundreds of ms).

### Spike-triggered average (STA)

The spatiotemporal spike-triggered average estimates the mean stimulus preceding spikes:

$$
\mathrm{STA}(\mathbf{x}, \tau) = \frac{1}{N_{\mathrm{sp}}} \sum_{k=1}^{N_{\mathrm{sp}}} S\big(\mathbf{x},\ t_k - \tau\big),
$$

where $t_k$ are spike times on the binned grid, $\tau \ge 0$ indexes time lag, and $N_{\mathrm{sp}}$ is the number of spikes used.

For spatial summaries (e.g., for center localization), the STA was reduced across time lags (e.g., by selecting or integrating over informative lags) and normalized to yield a spatial STA “strength” image.

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

## Static nonlinearity (LN output stage)

### Generator signal

Using the estimated linear spatiotemporal RF, we computed a scalar **generator signal** $g(t)$ by projecting the stimulus onto the filter:

$$
 g(t) = \sum_{\tau=0}^{T_h} \sum_{\mathbf{x}} K(\mathbf{x})\,h(\tau)\,S(\mathbf{x}, t-\tau),
$$

where $K(\mathbf{x})$ is a spatial RF (e.g., Gaussian fit or masked RF) and $h(\tau)$ is the temporal filter.

### Nonlinearity estimation

The static nonlinearity maps generator signal to firing rate. We estimated it empirically by binning $g(t)$ values and computing the mean firing rate within each bin:

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

Because responses were measured with repeated presentations, we quantified both:

- **Baseline reliability**: correlation between repeat-averaged responses across different repeat splits (or pairwise repeat correlations). This provides an empirical upper bound on explainable variance due to intrinsic trial-to-trial variability.

- **Prediction performance**: correlation between model prediction and measured response across time points, computed in a way that respects repeats (i.e., comparing predictions to the appropriate repeat-averaged response or using repeat-aware correlation measures).

Performance was summarized across stimulus conditions and then aggregated across cells/groups.

## Grouping and aggregation

Cells were grouped by response polarity (ON vs OFF; defined from the step-spot ON vs OFF phase response during patching) and by retinal location categories (e.g., nasal vs temporal). For each fitted quantity (spatial RF size/ellipticity, temporal filter parameters, nonlinearity parameters, center–surround indices, moving-bar model parameters and performance), we computed group means and uncertainty estimates (e.g., SEM across cells), and compared groups using appropriate statistical tests when needed.

## Practical reproducibility notes

- The pipeline relies on consistent stimulus–response alignment and careful spike detection; both materially affect STA estimates and downstream model fits.
- Regularization/constraints (e.g., fixed surround-to-center scale ratios for spot-size DoG fits; CSR penalties for moving-bar fits) are used to stabilize estimation given finite data.
- Wherever bootstrap resampling is used, confidence intervals should be reported from the bootstrap distribution (percentiles).
