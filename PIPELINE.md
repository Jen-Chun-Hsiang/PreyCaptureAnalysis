# PreyCaptureRGC analysis pipeline (reconstructed from code)

This repo is a MATLAB, script-driven analysis codebase for retinal ganglion cell (RGC) recordings during visual stimulation (moving noise/white noise, moving bars, spots/spot size).

There was no existing pipeline/README-style workflow doc found in the repo, so this document is a **code-grounded reconstruction** based on how scripts read/write data and on the explicit prerequisite checks embedded in the code.

> Note on style: many scripts are interactive/iterative (contain `keyboard;` and hard-coded recording lists). Treat this as a “what to run” map rather than a turnkey package.

---

## Folder conventions (as used by scripts)

### Inputs

- `Data\ephys\sections\...`  
  Per-recording ephys “section” `.mat` files (created by sectioning / matching scripts).

- `Data\Stimulation\Temporal_AlphaRGC_<recording>_<stimId>_Retina_1_<StimType>_1.mat`  
  Stimulus files. Commonly loaded variables include `OLED` (for pixel size, etc.).

- Excel metadata
  - `PatchRecordingAlphaRGCs_PreyCapture.xlsx` (used by white-noise aggregation scripts)
  - `PreyCaptureRGC_Data.xlsx` (used by `DataSummaryPlot.m`)

### Outputs

- `Results\MovingWhite\<recording>.mat`  
  Output of moving-noise/white-noise preprocessing per cell.

- `Results\MovingWhite\GaussianFitting_processed_*.mat`  
  Output bundle from `WhiteNoise_ONOFFalpha_Comparison.m` (spatial Gaussian RF fits, temporal filter fits, metadata, and often NL curves/params). Many downstream scripts **require** one specific processed bundle.

- `Results\MovingBar\*_moving_bar_processed.mat`  
  Output of moving bar preprocessing.

- `Results\MovingBar\*_moving_bar_LN_simulated_*.mat`  
  Output of LN simulation driven by the WN-fitted RF/TF parameters.

- `Results\MovingBar\*_moving_bar_fitted.mat`  
  Output of moving-bar fitting (correlations, fitted parameters).

- `Results\StationarySpot\...` and `Results\MovingSpot\...`  
  Outputs of stationary/moving spot preprocessing.

- `Results\Spots\SusAlpha_RawData.mat`  
  Input to the spot-size / center-surround analysis scripts in this repo (see below). This file is assumed to exist.

---

## Common numeric conventions

- Many pipelines downsample ephys/spikes to `Fz = 100` Hz for stimulus-response analyses.
- Time windows frequently used for stimulus averaging:
  - ON spot stimulus: `110:300`
  - OFF spot stimulus: `310:500`
  - (These are indices at 100 Hz; see spot-size scripts.)

---

## Pipeline overview (high level)

1. **(Optional/Upstream) Ephys↔Stim matching and ephys sectioning**
2. **Per-recording preprocessing** for each stimulus type (moving noise, moving bar, spots)
3. **White-noise aggregation / RF+TF+NL fitting** (central prerequisite)
4. **Downstream modeling & comparison** (moving bar LN simulation + fitting; spot size center/surround metrics)
5. **Summary/plot scripts**

Parallelization opportunities:
- Preprocessing steps are largely **per-recording independent** (safe to run in parallel on different recordings).
- Moving-bar LN simulation and fitting are also **per-recording independent**, once the WN processed bundle exists.

---

## Step 0 — Ephys↔Stim matching / sectioning (if needed)

Entry points:
- `DocumentEphys2StimulationData.m`

What it does (as used by other scripts):
- Produces/updates `EphysStimMatching.mat`.
- Can split recordings into “sections” (expects external `BasicFunctions` on MATLAB path).

Notes:
- Many preprocessing scripts assume the `Data\ephys\sections\...` files already exist.

---

## Step 1 — Moving noise / “white noise” preprocessing (per recording)

Entry point:
- `DocumentEphys2StimulationData_MovingNoise.m`

Inputs:
- One ephys section file
- One moving-noise stimulus file `Temporal_AlphaRGC_<rec>_<id>_..._MovingNoise_1.mat`

Outputs:
- `Results\MovingWhite\<recording>.mat` containing (typical):
  - `STAmat`, `stdSTA`, `masked_STAmat`, `masked_stdSTA`
  - Temporal RF estimate: `tRF`
  - Nonlinearity inputs: `PBs`, `FRs`

Notes:
- This script contains many per-recording spike detection thresholds and interactive logic.

---

## Step 2 — White-noise aggregation + RF/TF/NL fitting (central prerequisite)

Entry point:
- `WhiteNoise_ONOFFalpha_Comparison.m`

This is the core “spine” script: many downstream scripts explicitly error if its processed output is missing.

Inputs:
- Metadata spreadsheet: `PatchRecordingAlphaRGCs_PreyCapture.xlsx`
- All relevant per-recording files from `Results\MovingWhite\<recording>.mat`

Key computations:
- Spatial RF Gaussian fit (fits `stdSTA` → `gauss_est`)
- Temporal filter fit (difference-of-Gaussians via `GaussianTemporalFilter2.m` → `Gauss_TF_est`)
- Optional nonlinearity fitting (generator signal `PBs` vs firing `FRs`)

Outputs:
- `Results\MovingWhite\GaussianFitting_processed_*.mat`
  - Must contain, at minimum, `data_sets`, `cell_type`, `location` plus fitted parameter arrays.
  - Downstream moving-bar scripts commonly hardcode a particular filename (e.g. `GaussianFitting_processed_082025_1.mat`).

Related helper/plot scripts:
- `visualize_temporal_filter.m` (assumes variables from the main script are already in the workspace)
- `GetWhiteNoise_Nonlinearity.m` → produces `NL_params`, `NL_curves` (and figure outputs)
- `visualize_nonlinear_curve.m`
- `Nonlinearity_Analysis.m` (derives secondary NL metrics from `NL_params`)

---

## Step 3 — Spot size / center–surround analyses

There are two “layers” here:

### 3A) Spot-size summary/indices (group-level)

Entry points:
- `SpotSizeAnalysis_simple.m`
- `SpotSizeAnalysis_simple_byWavelength.m`

Inputs:
- `Results\Spots\SusAlpha_RawData.mat` which is expected to define a struct `a` with fields:
  - `x1` (time base)
  - `y_labels` (spot diameters)
  - group matrices like `DN_ONSus_RF_GRN`, etc.

Output(s):
- Figures under:
  - `Figures\SpotSizeSimple\`
  - `Figures\SpotSizeSimple_byWavelength\`

Computed metric (as implemented):
- **Size Index (SI)** per cell:
  - Choose “center” response `C` from the peak and its largest adjacent response among sizes `< 800`.
  - Choose “surround/large” response `S` as the mean response at diameters `[800, 1200]`.
  - `SI = (S - C) / C`.

Connection to moving-bar fitting:
- `MovingBar_Model_Prediction_Batch_Results.m` contains a note: “run `SpotSizeAnalysis_simple.m` to get the following values”, then hard-codes group-level constants like:
  - `ON_Nasal_CSR = -0.006`, `ON_Temporal_CSR = -0.09`, etc.
- These constants are later passed into moving-bar fitting as `CSR_value`.
  - Important sign convention: in `fit_linear_transform_with_surround.m`, the penalty uses `(gamma + CSR_value)^2` while enforcing `gamma >= 0`. So a **negative** `CSR_value` effectively targets `gamma ≈ -CSR_value`.

### 3B) DoG fitting of spot-size curves (per group and per cell)

Entry point:
- `SpotSizeAnalysis.m`

What it does:
- Computes mean response during stimulus per size, then:
  - fits a DoG model to the mean curve (either `fit_DoG_half_CDF` or `fitspotsizedog`)
  - optionally fits per-cell curves (spline interpolation + optional jitter)

Outputs:
- Figures under `Figures\SpotSizeFit\`
- Stores fit results in a MATLAB `results` struct in the workspace (the script primarily saves figures).

Helper:
- `fitspotsizedog.m` fits a constrained DoG-like form and reports a surround/center “weight ratio”:
  - `weight_ratio = (k_s*sigma_s^2) / (k_c*sigma_c^2)`

### 3C) RF center/surround parameters via bootstrap

Entry point:
- `RF_CenterSurround_Analysis.m`

What it does:
- Uses `estimate_rf_sizes.m` to estimate center/surround size parameters (bootstrapped) from the spot-size response curves.

Outputs:
- Saves into `Results\RF_Analysis\` (and creates figures).

---

## Step 4 — Moving bar pipeline (preprocess → LN simulate → fit)

### 4A) Preprocess moving bar recordings

Entry points:
- `DocumentEphys2StimulationData_MovingBar.m`
- `DocumentEphys2StimulationData_MovingBar_StageVSS.m` (variant that includes contrast dimension)

Outputs:
- `Results\MovingBar\<recording>_moving_bar_processed.mat` containing a `Data` array indexed by condition dimensions.

### 4B) Simulate LN response to moving bar using WN-derived filters

Entry points (several variants exist):
- `CompareModelPrediction8Responses_MovingBar_separateData.m` (and dated variants)
- `MovingBar_LinearNL_Simulation_steamline.m` (core simulator)

Prerequisite:
- A WN processed bundle must exist (e.g. `Results\MovingWhite\GaussianFitting_processed_082025_1.mat`).

Outputs:
- `Results\MovingBar\<recording>_<response>_moving_bar_LN_simulated_*.mat` containing `resp` (center drive) and `resp_s` (surround drive), plus the condition dimensions.

### 4C) Fit LN(+surround) transform and/or LNK variants

Entry points:
- `MovingBar_LinearNL_Fitting.m` (newer)
- `MovingBar_LinearNL_Fitting_082025.m` (older variant)

What it does:
- Loads processed moving bar data (`Data`) and matched LN simulation file (`resp`, `resp_s`).
- Builds a long trial vector and computes repeat reliability (`BaselineCorr`).
- Fits:
  - center-only linear transform: `fit_linear_transform.m`
  - center–surround transform: `fit_linear_transform_with_surround.m`
  - optional LNK variants through `test_LNK_fitting` (which uses functions like `fitLNK_rate_scw.m`).

CSR constraint support:
- If the variables `CSR_value` and `CSRStrength` exist in the workspace, `MovingBar_LinearNL_Fitting.m` passes them into `fit_linear_transform_with_surround`.
- Batch wrapper:
  - `MovingBar_Model_Prediction_Batch_Results.m` loads WN metadata and assigns a per-cell CSR constant based on ON/OFF × Nasal/Temporal.

Outputs:
- `Results\MovingBar\<recording>_moving_bar_fitted.mat` containing (typical):
  - `PredictionResults`, `BaselineCorr`
  - parameter structs like `LN_params_*`, `LNK_params_*`

### 4D) Summarize across cells

Entry point:
- `MovingBar_ONOFFalpha_Trace_Comparison.m`

What it does:
- Aggregates correlations/parameters across recordings.
- Cross-checks WN metadata alignment (cell type and location) against the WN processed bundle.

---

## Summary / exploratory scripts (often interactive)

- `DataSummaryPlot.m`
  - Loads many result types and computes additional indices.
  - Contains `keyboard;` breakpoints and is best treated as exploratory.

---

## “Which script should I run?” (recommended minimal path)

If your goal is the moving-bar modeling that depends on RF/TF estimates:

1. Preprocess moving noise per recording: `DocumentEphys2StimulationData_MovingNoise.m`
2. Aggregate and fit WN RF/TF (+NL optional): `WhiteNoise_ONOFFalpha_Comparison.m`
3. Preprocess moving bar per recording: `DocumentEphys2StimulationData_MovingBar.m` (or the `*_StageVSS.m` variant if your dataset has contrast)
4. Simulate LN for moving bar: `CompareModelPrediction8Responses_MovingBar_separateData.m` (or your preferred dated variant) → uses `MovingBar_LinearNL_Simulation_steamline.m`
5. Fit models per recording: `MovingBar_LinearNL_Fitting.m` (optionally via `MovingBar_Model_Prediction_Batch_Results.m`)
6. Summarize across cells: `MovingBar_ONOFFalpha_Trace_Comparison.m`

If your goal is spot-size center/surround quantification:

1. Ensure `Results\Spots\SusAlpha_RawData.mat` exists.
2. Run `SpotSizeAnalysis_simple.m` and/or `SpotSizeAnalysis.m`.
3. For bootstrapped center/surround estimates, run `RF_CenterSurround_Analysis.m`.

---

## Notes on legacy / multiple versions

There are many dated variants (e.g. `*_082025.m`, `*_092424.m`) that overlap in purpose.

Practical guidance:
- When scripts hardcode a processed WN bundle filename (e.g. `GaussianFitting_processed_082025_1.mat`), treat that as the “expected” output version for that downstream analysis.
- Prefer the non-dated entry points unless you know a given dataset requires the older logic.

---

## Common failure modes / troubleshooting

- “Run `WhiteNoise_ONOFFalpha_Comparison.m` first.”
  - Downstream moving-bar scripts will error if the WN processed bundle is missing.

- `keyboard;` stops execution
  - This is intentional debugging/interactive behavior in many scripts.

- Hard-coded network paths
  - Many scripts use absolute UNC paths; you may need to edit folder variables if your local mount differs.

---
