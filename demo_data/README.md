# Demo data

The `swr_cell_metrics_demo` package demonstrates the per-cell sharp-wave ripple metrics calculated by `compute_swr_cell_metrics`.

The package contains one real sleep-session LFP channel, selected sorted spike timestamp files, a precomputed `slow_cHFOs` event list, and the minimal timing and session metadata required by the analysis. It does not run SWR detection or include behavioural velocity data.

## Requirements

- MATLAB R2023.
- MATLAB Signal Processing Toolbox.
- MClust/Neuralynx readers that provide `readCRTsd`, `Data`, `Range`, `readSpikeDataOnly`, and `fixSpikes`.

Add the MClust location to `Matlab_code/classification/classification_config.json`, or pass it as `mclustPath` in the settings below.

## Run the demo

Run the following commands from the repository root. The commands create a separate runtime copy so the bundled input files remain unchanged.

```matlab
repoRoot = pwd;
sourceDemo = fullfile(repoRoot, 'demo_data', 'swr_cell_metrics_demo');
runtimeDemo = fullfile(repoRoot, 'demo_output', 'swr_cell_metrics_demo');

if exist(runtimeDemo, 'dir') == 7
    error('Choose a new runtime folder or remove the existing demo_output folder first.');
end

mkdir(fileparts(runtimeDemo));
copyfile(sourceDemo, fileparts(runtimeDemo));

loaded = load(fullfile(runtimeDemo, 'sessionInfo_template.mat'), 'sessInfo');
sessInfo = loaded.sessInfo;
sessInfo.mainDir = fullfile(runtimeDemo, 'session');
save(fullfile(runtimeDemo, 'sessionInfo.mat'), 'sessInfo');

settings = struct( ...
    'sessionInfoPath', fullfile(runtimeDemo, 'sessionInfo.mat'), ...
    'allCellsPath', fullfile(runtimeDemo, 'All_Cells_combined.mat'), ...
    'sleepFolders', {{'s1'}}, ...
    'mclustPath', '/path/to/MClust');

result = compute_swr_cell_metrics(settings);
```

Use the `mclustPath` value only when the MClust path is not already set in `classification_config.json`.

## Expected output

The command reports one processed session and creates:

```text
demo_output/swr_cell_metrics_demo/All_Cells_combined.mat
demo_output/swr_cell_metrics_demo/session/s1/processedData/SWR_data.mat
```

`All_Cells_combined.mat` contains the per-cell fields `S1_PSP`, `S1_SFI`, and `S1_SpPR`, as well as additional S1 SWR metrics. `SWR_data.mat` contains the event-by-cell matrices and intermediate quantities used to calculate them.

The packaged `manifest.txt` records the number of input SWR events and included units. Check that `result.sessions.status` is `processed` and that `result.sessions.swr_counts(1)` matches the manifest event count.

Measure and record the demo runtime on a standard desktop before the repository release.
