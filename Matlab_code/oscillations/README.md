# Oscillations

This folder contains the SWR/HFO detection step that generates the per-sleep event files used by the cell-level SWR metrics.

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
