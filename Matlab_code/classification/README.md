# Cell Classification Workflow

This folder contains the code that produces the cell-type labels used downstream.

Pipeline order:

1. `CCG_based_classifier_modified.m`
   - Generates CCG review-session files for excitatory/inhibitory monosynaptic candidate pairs.
   - Uses `classification_config.json` for local CellExplorer and MClust dependency paths.
   - Review/curation outputs are used as the ground-truth source for the next step.

2. `review_CCG_connection_labels.m`
   - Opens the manual CCG review UI for accepting, discarding, or relabeling candidate connections.
   - Writes curated session files and accepted-connection summaries.

3. `rebuild_All_Cells_combined_CCG_classification_from_curated.m`
   - Rebuilds `All_Cells_combined.CCGbased_CellClassfication` from curated CCG review files.
   - CCG labels use `1` for excitatory/principal ground truth and `2` for inhibitory/interneuron ground truth.

4. `GMM_based_classifications.m`
   - Trains and applies the principal/interneuron GMM using CCG-based ground truth labels.
   - The current feature set is spatial coverage, ACG mean, and classification firing rate.
   - Writes `All_Cells_combined.GMM_based_classification_days`.

5. `Classifications_Spatial_Info.m`
   - Combines `GMM_based_classification_days` with place-field outputs.
   - Writes `final_classification_numeric`:
     - `1`: principal non-place cell
     - `2`: place cell
     - `3`: interneuron
     - `4`: unclassified

Use `classification_config.json` to store local dependency paths when running code that needs CellExplorer or MClust.
