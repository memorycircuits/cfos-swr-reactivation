# Oscillations

This folder contains LFP-based oscillation analyses used by later cell-level metrics and exports.

## SWR/HFO detection

Run:

```matlab
detect_swr_events()
```

By default this processes `s1` and `s2` for all entries in `sessInfo` and writes:

```text
<session>/s1/processedData/_allE_numSD3.5_HighPwrCycles4.mat
<session>/s2/processedData/_allE_numSD3.5_HighPwrCycles4.mat
```

The default threshold is the settled SWR detection setting:

```matlab
numSD = 3.5;
numHighPwrCycles = 4;
peakSD = 3;
```

`detect_swr_events.m` is the cleaned replacement for `total_swr_hj.m`. `detect_swr_events_blanco.m` contains the Blanco-style RMS and peak-threshold event detector that was previously in `blanco_hj.m`.

The generated `.mat` files contain `slow_cHFOs`, `fast_cHFOs`, `all_cHFOs`, `lfp_bp`, `HFO_rate`, quiet/fast time summaries, and the `detection_threshold` structure.

Dependencies:

- `sessInfo`, usually from `sessionInfo.mat`.
- MClust/Neuralynx readers for `readCRTsd`, `Data`, and `Range`.
- The `detectHFOs` package for `blanco_bp` and `start_stop`.
- `findpeaksmine`.
- MATLAB Signal Processing Toolbox functions used by the original detector (`filtfilt`, `findpeaks`, `downsample`).

Store the MClust path in `../classification/classification_config.json`. Extra dependency folders can be passed without editing the file:

```matlab
detect_swr_events(struct( ...
    'oscillationAnalysisPath', '/path/to/detectHFOs', ...
    'additionalPaths', {{'/path/to/findpeaksmine'}}));
```

Useful options:

```matlab
detect_swr_events(struct( ...
    'sessionInfoPath', '/path/to/sessionInfo.mat', ...
    'sessionIndices', 1:10, ...
    'sleepFolders', {{'s1', 's2'}}, ...
    'overwriteExisting', false));
```

## SWR phase analysis

Run after SWR/HFO detection:

```matlab
compute_swr_phase_analysis()
```

This is the cleaned SWR spike-phase analysis from `instantaneous_phase_v1.m`. `SWR_continuous_rate_mh.m` was checked and is a separate continuous SWR-rate helper, not the phase analysis.

The phase analysis reads each sleep session's SWR event file:

```text
<session>/s1/processedData/_allE_numSD3.5_HighPwrCycles4.mat
<session>/s2/processedData/_allE_numSD3.5_HighPwrCycles4.mat
```

It writes per-session phase details to:

```text
<session>/s1/processedData/Instantaneous_phase.mat
<session>/s2/processedData/Instantaneous_phase.mat
```

It also updates `All_Cells_combined` with the downstream fields:

```text
S1_cells_phase_mean_angle
S1_cells_phase_R
S1_cells_phase_p_val
S1_cells_phase_z
S1_cells_phase_mean_angle_troughAligned
S1_cells_phase_R_troughAligned
S1_cells_phase_p_val_troughAligned
S1_cells_phase_z_troughAligned
```

The same fields are written with the `S2_` prefix for `s2`.

Default settings:

```matlab
swrFileName = '_allE_numSD3.5_HighPwrCycles4.mat';
phaseFileName = 'Instantaneous_phase.mat';
phaseDownsamplingFactor = 4;
sleepFolders = {'s1', 's2'};
```

Dependencies:

- SWR event files from `detect_swr_events.m`.
- `All_Cells_combined`.
- MClust/Neuralynx readers for `readCRTsd`, `Data`, and `Range`.
- Spike-loading helpers for `readSpikeDataOnly` and `fixSpikes`.
- The `detectHFOs` package for `blanco_bp`.
- MATLAB Signal Processing Toolbox functions (`hilbert`, `filtfilt`, `findpeaks`, `downsample`).

## Open-field theta phase encoding

Run:

```matlab
compute_open_field_theta_phase_encoding()
```

This is the cleaned public-repo version of `theta_phase_encoding_of.m`. It assigns instantaneous theta-band LFP phase to open-field spikes during running periods and writes per-cell circular phase metrics:

```text
All_Cells_combined(iii).thetaPhase_pref
All_Cells_combined(iii).thetaPhase_R
All_Cells_combined(iii).thetaPhase_p
All_Cells_combined(iii).thetaPhase_z
All_Cells_combined(iii).thetaPhase_nSpikes
```

It also saves per-session details to:

```text
<session>/processedData/OF_theta_phase.mat
```

Default theta settings:

```matlab
thetaBand = [6 14];
speedThreshold = 2;
minSpikes = 5;
openFieldFolders = {'of1', 'of2', 'of3'};
```

Example with explicit paths:

```matlab
compute_open_field_theta_phase_encoding( ...
    'sessionInfoPath', '/path/to/sessionInfo.mat', ...
    'allCellsPath', '/path/to/All_Cells_combined.mat', ...
    'additionalPaths', {'/path/to/spike-reader-helpers'});
```

Dependencies:

- `sessInfo`, usually from `sessionInfo.mat`.
- `All_Cells_combined`.
- MClust/Neuralynx readers for `readCRTsd`, `Data`, and `Range`.
- Spike-loading helpers for `readSpikeDataOnly` and `fixSpikes`.
- MATLAB Signal Processing Toolbox functions (`hilbert`, `butter`, `filtfilt`, `downsample`).
