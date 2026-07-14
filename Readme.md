# cFos SWR Reactivation
## Overview

This repository contains code associated with the study "Inhibitory and excitatory cFos engram neurons are preferentially reactivated by sharp wave ripples" by Javed, Robles-Hernandez et al. from the Memory Circuits Lab. Here we investigated how cFos-tagged hippocampal CA1 neurons are recruited during sharp wave ripples.

In this project, we combined cFos-dependent tagging, chronic in vivo electrophysiology, local field potential recordings, optotagging, spatial coding analyses, interneuron characterization, and sharp wave ripple analyses to ask whether experience-associated cFos-tagged neurons are preferentially recruited during post-behavior rest.

The study focuses on both excitatory and inhibitory components of a cFos-tagged spatial ensemble, including cFos-tagged place cells, non-place principal cells, and interneurons.

## System requirements

The source code is written in MATLAB. It was tested with a MATLAB 2023 release.

- **MATLAB:** MATLAB R2023 (tested).
- **MATLAB toolboxes:** Signal Processing Toolbox and Statistics and Machine Learning Toolbox.
- **Third-party MATLAB dependencies:** MClust/Neuralynx readers and CellExplorer (for CCG classification).
- **Operating system:** The code is intended for standard desktop operating systems supported by MATLAB.
- **Hardware:** No non-standard hardware is required to run the analyses. Memory and storage requirements depend on the size of the electrophysiology recordings.

## Installation

1. Download or clone this repository.
2. Install MATLAB R2023 and the required toolboxes.
3. Install the third-party dependencies listed above and record their local paths in [`Matlab_code/classification/classification_config.json`](Matlab_code/classification/classification_config.json). At minimum, configure `cellExplorerPath` and `mclustPath`; pass the `detectHFOs` and `findpeaksmine` folders through the `additionalPaths` setting when running oscillation analyses.
4. Start MATLAB and add the repository code to the MATLAB path:

```matlab
addpath(genpath('/path/to/cfos-swr-reactivation/Matlab_code'))
```

The MATLAB-path setup takes less than one minute. Once MATLAB and the third-party dependencies are already installed, configuration normally takes less than 10 minutes; installing external dependencies may take longer.

## Demo data

The [`demo_data`](demo_data/README.md) folder documents how to run the packaged SWR cell-metrics example. It calculates per-cell SWR participation, firing-rate increase, and spikes-per-ripple fields from a small real sleep-session subset with precomputed SWR events.

## Running the analyses

The electrophysiology analyses use a `sessionInfo.mat` file containing `sessInfo` and an `All_Cells_combined.mat` file. Pass both locations explicitly when they are not in the repository root or `Data` folder:

```matlab
settings = struct( ...
    'sessionInfoPath', '/path/to/sessionInfo.mat', ...
    'allCellsPath', '/path/to/All_Cells_combined.mat');
```

### Sharp-wave ripple and phase analyses

Run SWR detection first, followed by the desired cell-level analyses:

```matlab
detect_swr_events(settings)
compute_swr_cell_metrics(settings)
compute_swr_phase_analysis(settings)
```

For open-field theta-phase encoding, run:

```matlab
compute_open_field_theta_phase_encoding( ...
    'sessionInfoPath', settings.sessionInfoPath, ...
    'allCellsPath', settings.allCellsPath)
```

See [`Matlab_code/oscillations/README.md`](Matlab_code/oscillations/README.md) for required session files, output filenames, and optional settings.

### Cell classification and place-cell metrics

Compute the GMM features before running the GMM classifier:

```matlab
compute_gmm_cell_properties(settings)
GMM_based_classifications(settings)
Classifications_Spatial_Info(settings)
compute_place_cell_reuse_metrics(settings)
```

The GMM classifier requires curated CCG ground-truth labels. To generate and review these labels, run `CCG_based_classifier_modified(settings)`, curate the generated review sessions with `review_CCG_connection_labels`, and then run `rebuild_All_Cells_combined_CCG_classification_from_curated`. The complete sequence is described in [`Matlab_code/classification/README.md`](Matlab_code/classification/README.md).

### Optotagging and CSV exports

Run the SALT optotagging analysis before exporting its summaries:

```matlab
SALT_opto_analysis(settings)
CSV_file_export_first(settings)
```

### Image quantification

Provide the path to the quantification workbook and an output folder:

```matlab
plot_CB_quantification_figures( ...
    '/path/to/Quantifications.xlsx', ...
    'OutputFolder', '/path/to/output/CB')

plot_interneuron_quantification_figures( ...
    '/path/to/Quantifications.xlsx', ...
    'OutputFolder', '/path/to/output/interneurons')
```

The figure functions write figures and summary tables to their specified output folders. The cell-property workflow and its input expectations are documented in [`Matlab_code/cell_properties/README.md`](Matlab_code/cell_properties/README.md).

## Preprint

**bioRxiv**: https://www.biorxiv.org/content/10.1101/2024.12.17.628897


## Lab

This work was carried out in the Memory Circuits Lab. <br> <br>
<img src="lab/memorycircuits_lab_logo.png" alt="Memory Circuits Lab logo" width="250"/> <br>

At the German Center for Neurodegenerative Diseases  (DZNE) <br>
and  <br>
Charité - Universitätsmedizin Berlin <br>
Berlin, Germany <br>

## Contact

For questions about the code or analyses, please contact: <br> <br>

Silvia Viana da Silva <br>
Memory Circuits Lab <br>
German Center for Neurodegenerative Diseases (DZNE) <br>
Berlin, Germany <br> <br>

or  <br>
Matthias Haberl <br>
Charité - Universitätsmedizin Berlin <br>
Berlin, Germany <br>


## License

Please see the repository license for terms of use. <br>
