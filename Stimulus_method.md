# Visual stimulus methods (Cogent + Stage/LightCrafter)

This document summarizes how visual stimuli were generated and presented for ex vivo retina electrophysiology. It is written for inclusion in a peer‑review paper Methods section: **what the stimuli were, how they were parameterized, and how timing/contrast were controlled**, independent of any downstream analysis scripts.

The stimulus code lives primarily in the [Running_scripts/](Running_scripts/) folder (experiment “launcher” scripts) and supporting functions in the repository root and [LightCrafterProjector/](LightCrafterProjector/).

## General stimulus timing and synchronization

### Frame timing

Two presentation backends are used:

1. **Cogent Graphics (MATLAB)** for OLED-based stimuli (checkerboard noise, spots, moving spots). Timing is driven by `cgflip`, with frame timestamps recorded.
2. **Stage Visual Stimulation System (Stage VSS)** for LightCrafter-based moving bar stimuli. Timing is driven by the Stage `Presentation` engine, and per-presentation `flip_times` are retrieved from the Stage client.

### Synchronization to electrophysiology

Stimulus onset/offset markers are emitted as TTL pulses on a digital output line (parallel port / DIO), using functions such as `openparallel`, `pin16high`, `pin16low`, and (for some protocols) `UniqueTrigger`. These pulses provide a hardware synchronization signal to align stimulus epochs to recorded electrophysiology.

### Intensity calibration / gamma correction

For OLED/Cogent spot and moving‑spot stimuli, commanded intensities were converted to device pixel values using an empirically measured intensity calibration curve stored in `NF_IntensityMeasurement_032921_NF_2_002.mat`.

The mapping is implemented by `ProjIntCalib.m`, which inverts a fitted sigmoidal intensity function to compute the “computer intensity” that produces a desired target intensity/contrast on the display.

## Display geometry and units

Stimulus parameters are specified in **microns on the retina** and converted to screen pixels using a calibrated `pixelSize` (µm/pixel). The calibrated `pixelSize` is stored alongside each stimulus file (e.g., as `OLED.pixelSize` or `PROJ.pixelSize`).

Typical values used in scripts include:

- OLED canvas: `width = 800` px, `height = 600` px, refresh `60` Hz.
- Calibrated pixel size on retina depends on optical configuration; scripts include values such as `2.0`, `2.5`, and `4.375` µm/pixel.

For Stage/LightCrafter protocols, a `center_position` in projector pixels is defined to align stimulus coordinates to the retinal field.

---

## Stimulus 1: Moving noise (jittered checkerboard “white noise”)

**Goal:** estimate fine spatiotemporal receptive field structure while using relatively large checkerboard elements that reliably drive spikes.

### Spatial pattern

At each stimulus update, the displayed image is a binary checkerboard of independent samples:

- A stimulus grid of size `nRows × nCols` is computed from the display size and the requested checker size.
- Each grid element is sampled independently as `0` or `1` (black/white) using a seeded pseudo-random generator.

If the checker size is `NoiseGridSize` (µm) and pixel size is `pixelSize` (µm/pixel), then:

$$
 nRows = \left\lceil \frac{H_{px}}{\mathrm{NoiseGridSize}/\mathrm{pixelSize}} \right\rceil,\qquad
 nCols = \left\lceil \frac{W_{px}}{\mathrm{NoiseGridSize}/\mathrm{pixelSize}} \right\rceil.
$$

### Temporal update rate and duration

The checkerboard is updated at a fixed rate `NoiseFz` (Hz) for a total duration `minT` (minutes). The presentation is split into blocks of length `NoiseBlockLength` (seconds):

$$
 nBlock = \left\lceil \frac{60\,\mathrm{minT}}{\mathrm{NoiseBlockLength}} \right\rceil,\qquad
 nFrame = \lceil \mathrm{NoiseFz}\cdot\mathrm{NoiseBlockLength} \rceil.
$$

Typical parameters (example from launch scripts):

- `NoiseGridSize = 100 µm`
- `NoiseFz = 10 Hz`
- `NoiseBlockLength = 10 s`
- `minT ≈ 10–12 min`

### Jitter (“moving noise”)

To increase effective spatial sampling, the **entire checkerboard is randomly jittered each update**. For each frame, independent horizontal/vertical offsets are sampled from a discrete set `MovingRange` (µm) and converted to pixels:

$$
 \Delta x_t,\Delta y_t \sim \mathrm{UniformDiscrete}(\mathrm{MovingRange}/\mathrm{pixelSize}).
$$

The checkerboard sprite is drawn at `(Δx_t, Δy_t)` relative to the screen origin.

### Repeated vs. unique blocks (frozen noise segments)

A fraction `RepeatProb` of blocks are designated as repeated (“frozen”) blocks to allow repeat‑based reliability estimation.

- **Unique blocks** are generated from a pseudo‑random stream seeded by `NoiseUniqueSeed` and advanced continuously across frames.
- **Repeated blocks** are generated from a stream seeded by `NoiseRepeatSeed` and *restarted* so the same frame sequence is replayed.

Typical example parameters:

- `RepeatProb = 0.125`
- `NoiseUniqueSeed = 119`
- `NoiseRepeatSeed = 78`

### Logged metadata

The moving-noise code records a frame table including frame flip timestamps, block index, a repeated/unique flag, and frame number. Jitter offsets (`MovingSteps`) are also stored.

### Optional dual-color / dual-region moving noise

A dual-sprite variant exists (`MovingWhiteNoise_SciRig_ephys_UVGreen.m`) that loads two checkerboard sprites into different color channels and draws them at two configurable positions (`ColorPos = [x1 y1 x2 y2]`). One channel can be attenuated by a scalar factor (e.g., `GreenAttenuation`). This supports simultaneous stimulation of two spatial regions and/or spectral channels.

---

## Stimulus 2: Stationary spots (size × contrast; ON and OFF steps)

Stationary spot stimuli were used both for center–surround characterization and for determining ON vs OFF polarity.

### Geometry

- Spots are circular disks centered at a specified location (`Center`).
- Spot diameters are specified in µm and converted to pixel radius:

$$
 r_{px} = \frac{d/2}{\mathrm{pixelSize}}.
$$

A circular **aperture/mask** can be applied (e.g., `MaskSize` µm in `Spot_SciRig_ephys_051724.m`) by drawing a background disk and then drawing the spot within it.

### Contrast/intensity specification

Two contrast conventions appear in the scripts:

1. **Absolute intensity list** (`StationalSpot_SciRig_ephys.m`): `IN.contrasts` is treated as a set of target intensities (typically spanning 0–1), mapped through `ProjIntCalib`.

2. **Increment/decrement around mean gray** (`Spot_SciRig_ephys_051724.m`): contrasts are specified as amplitudes around a mean of 0.5.

For contrast amplitude $c$:

$$
 I_{ON} = 0.5 + \frac{c}{2},\qquad I_{OFF} = 0.5 - \frac{c}{2}.
$$

These target values are then converted to device pixel commands via `ProjIntCalib`.

### Timing

Two common timing structures are used:

- **Single-phase spot with ISI** (`StationalSpot_SciRig_ephys.m`):
  - Present spot for `presentT` seconds.
  - Return to mean background for `ISI` seconds.

- **Explicit ON and OFF phases** (`Spot_SciRig_ephys_051724.m`):
  - Present ON spot for `presentT(1)` seconds.
  - Present OFF spot for `presentT(2)` seconds.
  - Inter-stimulus interval `ISI` between size/contrast conditions.

### Randomization and repeats

Size×contrast conditions are enumerated and randomized per global repeat. Scripts support both “local” repeats within each condition and “global” repeats over the full condition set.

### Logged metadata

The spot scripts store a table containing stimulus identifiers (diameter index, contrast index, ON vs OFF phase) and the flip timestamp for each transition.

---

## Stimulus 3: Moving spots (diameter × speed × direction)

Moving spots are circular disks translated across the retinal field.

### Geometry and motion

- Spot diameter `d` (µm) determines radius $r_{px}$ as above.
- Motion is along directions given in degrees (`directions`, e.g., 0°, 90°, 180°).
- Each epoch starts with the spot positioned at a distance `startposition` (µm) from the center, then moved across the center and beyond for a total travel distance of approximately `2*(startposition + radius)`.

### Speeds

Speeds are specified in µm/s (e.g., 500–8000 µm/s) and converted to pixels/s:

$$
 v_{px/s} = \frac{v_{\mu m/s}}{\mathrm{pixelSize}}.
$$

### Timing and repeats

For each diameter×speed×direction condition:

- The spot is rendered continuously during motion.
- An ISI (`ISI`, seconds) is inserted between trials.
- Conditions are randomized within each repeat.

### Intensities

Spot and background intensities are specified by `spotcontrast` and `bgcontrast`, mapped through `ProjIntCalib`.

---

## Stimulus 4: Moving bars (Stage VSS + LightCrafter)

Moving bars were presented using the **Stage Visual Stimulation System** driving a **Texas Instruments LightCrafter 4500** in pattern mode.

### Hardware configuration and low-bit / high-frame-rate presentation

Launcher scripts (e.g., [Running_scripts/MovingBar_StageVSS_Script100224.m](Running_scripts/MovingBar_StageVSS_Script100224.m)) configure:

- `PROJ.pixelSize` (µm/pixel on retina; e.g., 4.375 µm/pixel)
- LED channel selection (`LED_color_str`, e.g., blue or green)
- Maximum LED current (`max_led_current`)
- A *pattern up-factor* (`target_frame_up_factor`, e.g., 12) to increase effective pattern rate.

Stage is configured with a `PatternRenderer` and the LightCrafter is set to pattern mode with attributes chosen to match the desired pattern rate.

### Stimulus parameterization

Each moving-bar trial is parameterized by:

- **Bar width** $w$ (µm): e.g., `[50, 100, 200, 400, 800]`
- **Bar speed** $v$ (µm/s): e.g., `[500, 1000, 2000, 4000, 8000]`
- **Direction** $\theta$ (deg): e.g., `[0, 90, 180]`
- **Bar height** (µm): e.g., `Height = 800` µm (rendered as a long rectangle)
- **Start distance** (µm): `startposition` (how far from center the bar begins)
- **Background disc**: a circular region of radius `BGdiscradius` (µm) centered on the stimulus center
- **Background intensity set**: `BGcontrasts` (e.g., `[0, 0.33, 0.66]`)
- **Bar polarity**: black or white bars (`BW = 'B'` or `'W'`)

All distances are converted to pixels using `PROJ.pixelSize`.

### Geometry and motion

For each trial:

1. The bar is initialized at a position offset from the stimulus center such that the bar starts outside the center by `startposition + w/2`.
2. The bar is translated at constant velocity along the specified direction.
3. The total stimulus duration is set so the bar traverses from one side to the other:

$$
 T_{stim} = \frac{2\,(startposition + w/2)}{v}.
$$

A Stage `PropertyController` updates the bar position as a function of time, ensuring constant-velocity motion.

### Contrast / polarity conventions

In the “black bar” convention (`BW='B'`):

- Bar intensity is set to the minimum (black).
- Background intensity is chosen from `1 - BGcontrasts` so that larger `BGcontrasts` correspond to darker backgrounds.

In the “white bar” convention (`BW='W'`):

- Bar intensity is set to the maximum (white).
- Background intensity is chosen directly from `BGcontrasts`.

A circular background disc is rendered at the specified background intensity; the Stage canvas background outside the disc is set to 0.

### Trial structure, randomization, and repeats

- All combinations of width×speed×direction are enumerated and randomized within each repeat.
- Background intensities are also randomized within each repeat.
- Each trial is followed by an inter-stimulus interval `ISI` (seconds).

### Logged metadata

For each trial, the code records:

- repeat index
- event type (start/end/ISI markers)
- width/speed/direction indices
- background contrast index
- Stage-reported flip timestamps (`flip_times(1)` and `flip_times(end)`).

---

## Notes on script versions and missing dependencies

The `Running_scripts/` folder contains multiple dated versions of launchers (e.g., `*_090524.m`, `*_100224.m`). For documentation purposes, the **latest versions** were treated as canonical for parameter definitions:

- Stage moving bars: `StageVSS_MovingBar_100224.m` + `MovingBar_StageVSS_Script100224.m`.

Some older OLED “Green” wrappers reference helper functions that are not present in this repository (e.g., `MovingWhiteNoise_SciRig_ephys_Green`, `Spot_SciRig_ephys_Green`). Where applicable, stimulus behavior is described using the closest available implementations (`MovingWhiteNoise_SciRig_ephys.m`, `MovingWhiteNoise_SciRig_ephys_UVGreen.m`, `StationalSpot_SciRig_ephys.m`, `Spot_SciRig_ephys_051724.m`).
