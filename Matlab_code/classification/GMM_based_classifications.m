function GMM_based_classifications(custom_settings)
% Train and apply principal/interneuron GMM labels from CCG-curated ground truth.

if nargin < 1 || isempty(custom_settings)
    custom_settings = struct();
end

all_cells_path = char(string(get_override_value(custom_settings, 'allCellsPath', '')));
if isempty(all_cells_path)
    all_cells_path = resolve_all_cells_path();
end

output_root = char(string(get_override_value(custom_settings, 'outputRoot', fullfile(pwd, 'CellClassification'))));
figure_outfolder = char(string(get_override_value(custom_settings, 'figureOutfolder', fullfile(output_root, 'figures'))));
curve_outfolder = char(string(get_override_value(custom_settings, 'curveOutfolder', fullfile(output_root, 'classification_curves'))));
save_updated_all_cells = logical(get_override_value(custom_settings, 'saveUpdatedAllCells', true));

if ~exist(figure_outfolder, 'dir')
    mkdir(figure_outfolder);
end
if ~exist(curve_outfolder, 'dir')
    mkdir(curve_outfolder);
end

loaded = load(all_cells_path, 'All_Cells_combined');
if ~isfield(loaded, 'All_Cells_combined')
    error('All_Cells_combined was not found in %s.', all_cells_path);
end
All_Cells_combined = loaded.All_Cells_combined;

[classification_table, feature_matrix, feature_names, ccg_labels] = collect_gmm_training_table(All_Cells_combined);
if isempty(feature_matrix)
    error('No cells had the complete feature set needed for GMM classification.');
end
if ~any(ccg_labels == 1) || ~any(ccg_labels == 2)
    error('CCG ground-truth labels must include both excitatory/principal and inhibitory/interneuron cells.');
end

write_single_feature_probability_curves(classification_table, feature_names, curve_outfolder);

GMM_based_classification = make_classification_GMM(feature_matrix, feature_names, ccg_labels, figure_outfolder);
All_Cells_combined = assign_gmm_classification_to_sessions(All_Cells_combined, classification_table, GMM_based_classification);

if save_updated_all_cells
    save(all_cells_path, 'All_Cells_combined', '-v7.3');
end

save(fullfile(output_root, 'GMM_cell_classification_results.mat'), ...
    'classification_table', 'feature_matrix', 'feature_names', 'ccg_labels', 'GMM_based_classification');

fprintf('Saved GMM principal/interneuron classification results to %s\n', output_root);
fprintf('Feature set: %s\n', strjoin(feature_names, ', '));
fprintf('CCG-labeled training cells: %d\n', nnz(ccg_labels == 1 | ccg_labels == 2));
fprintf('Cells classified with complete features: %d\n', numel(GMM_based_classification));
end

function [classification_table, X_all, feature_names, ccg_labels] = collect_gmm_training_table(All_Cells_combined)

feature_names = {'Spatial coverage', 'ACG mean', 'Firing rate'};
SessionIndex = [];
CellIndex = [];
SpatialCoverage = [];
ACGMean = [];
FiringRate = [];
CCGLabel = [];

for session_index = 1:numel(All_Cells_combined)
    session_entry = All_Cells_combined(session_index);
    n_cells = infer_session_unit_count(session_entry);
    if n_cells == 0
        continue
    end

    spatial_coverage = get_numeric_field(session_entry, 'spatial_coverage_meanOFs', n_cells, NaN);
    acg_mean = get_numeric_field(session_entry, 'acg_mean', n_cells, NaN);
    firing_rate = get_numeric_field(session_entry, 'classific_firingRate', n_cells, NaN);
    ccg_label = get_numeric_field(session_entry, 'CCGbased_CellClassfication', n_cells, 0);

    complete_feature_mask = isfinite(spatial_coverage) & isfinite(acg_mean) & isfinite(firing_rate);
    if ~any(complete_feature_mask)
        continue
    end

    cell_indices = find(complete_feature_mask);
    SessionIndex = [SessionIndex; repmat(session_index, numel(cell_indices), 1)]; %#ok<AGROW>
    CellIndex = [CellIndex; cell_indices(:)]; %#ok<AGROW>
    SpatialCoverage = [SpatialCoverage; spatial_coverage(cell_indices)]; %#ok<AGROW>
    ACGMean = [ACGMean; acg_mean(cell_indices)]; %#ok<AGROW>
    FiringRate = [FiringRate; firing_rate(cell_indices)]; %#ok<AGROW>
    CCGLabel = [CCGLabel; ccg_label(cell_indices)]; %#ok<AGROW>
end

classification_table = table(SessionIndex, CellIndex, SpatialCoverage, ACGMean, FiringRate, CCGLabel);
X_all = [SpatialCoverage(:), ACGMean(:), FiringRate(:)];
ccg_labels = CCGLabel(:);
end

function write_single_feature_probability_curves(classification_table, feature_names, curve_outfolder)

feature_columns = {'SpatialCoverage', 'ACGMean', 'FiringRate'};

for iFeature = 1:numel(feature_columns)
    prop_name = feature_columns{iFeature};
    prop_values = classification_table.(prop_name);
    class_labels = classification_table.CCGLabel;

    valid = isfinite(prop_values) & (class_labels == 1 | class_labels == 2);
    prop_values = prop_values(valid);
    class_labels = class_labels(valid);

    data_exc = prop_values(class_labels == 1);
    data_inh = prop_values(class_labels == 2);
    if isempty(data_exc) || isempty(data_inh)
        warning('Skipping %s because one class has no valid data.', prop_name);
        continue
    end

    prior_exc = numel(data_exc) / (numel(data_exc) + numel(data_inh));
    prior_inh = numel(data_inh) / (numel(data_exc) + numel(data_inh));

    x_min = min(prop_values);
    x_max = max(prop_values);
    x_range = max(x_max - x_min, 1);
    lower = x_min - 0.08 * x_range;
    upper = x_max + 0.08 * x_range;
    support_arg = 'unbounded';

    if strcmp(prop_name, 'SpatialCoverage')
        lower = 0;
        support_arg = [0 x_max];
    elseif strcmp(prop_name, 'FiringRate') || strcmp(prop_name, 'ACGMean')
        lower = max(0, lower);
        support_arg = 'positive';
    end

    x_fine = linspace(lower, upper, 500);
    [f_exc, x_eval] = ksdensity(data_exc, x_fine, 'Support', support_arg, 'BoundaryCorrection', 'reflection');
    f_inh = ksdensity(data_inh, x_fine, 'Support', support_arg, 'BoundaryCorrection', 'reflection');

    p_exc_given_x = (f_exc .* prior_exc) ./ ((f_exc .* prior_exc) + (f_inh .* prior_inh));
    p_inh_given_x = (f_inh .* prior_inh) ./ ((f_exc .* prior_exc) + (f_inh .* prior_inh));

    p_exc_given_x(~isfinite(p_exc_given_x)) = NaN;
    p_inh_given_x(~isfinite(p_inh_given_x)) = NaN;

    fig = figure('Color', 'w', 'Position', [100 100 900 700]);
    ax = axes('Parent', fig);
    hold(ax, 'on');
    plot(ax, x_eval, p_exc_given_x, '-', 'Color', [0.2 0.4 0.8], 'LineWidth', 2.5);
    plot(ax, x_eval, p_inh_given_x, '-', 'Color', [0.8 0.2 0.2], 'LineWidth', 2.5);
    scatter(ax, data_exc, 0.02 * ones(size(data_exc)), 16, [0.2 0.4 0.8], 'filled', 'MarkerFaceAlpha', 0.35);
    scatter(ax, data_inh, 0.98 * ones(size(data_inh)), 16, [0.8 0.2 0.2], 'filled', 'MarkerFaceAlpha', 0.35);
    xlabel(ax, feature_names{iFeature}, 'Interpreter', 'none');
    ylabel(ax, 'P(class | value)');
    title(ax, sprintf('CCG ground-truth class probability: %s', feature_names{iFeature}), 'Interpreter', 'none');
    legend(ax, {'Principal/excitatory', 'Interneuron/inhibitory'}, 'Location', 'best', 'Box', 'off');
    xlim(ax, [lower upper]);
    ylim(ax, [0 1]);
    grid(ax, 'on');
    exportgraphics(fig, fullfile(curve_outfolder, sprintf('Distribution_%s.pdf', prop_name)), ...
        'ContentType', 'vector', 'BackgroundColor', 'white');
    close(fig);

    writematrix(x_eval(:), fullfile(curve_outfolder, sprintf('classification_%s_x_fine.csv', prop_name)));
    writematrix(p_exc_given_x(:), fullfile(curve_outfolder, sprintf('classification_%s_p_principal.csv', prop_name)));
    writematrix(p_inh_given_x(:), fullfile(curve_outfolder, sprintf('classification_%s_p_interneuron.csv', prop_name)));
end
end

function All_Cells_combined = assign_gmm_classification_to_sessions(All_Cells_combined, classification_table, GMM_based_classification)

for session_index = 1:numel(All_Cells_combined)
    n_cells = infer_session_unit_count(All_Cells_combined(session_index));
    All_Cells_combined(session_index).GMM_based_classification_days = zeros(n_cells, 1);
end

for iRow = 1:height(classification_table)
    session_index = classification_table.SessionIndex(iRow);
    cell_index = classification_table.CellIndex(iRow);
    if session_index < 1 || session_index > numel(All_Cells_combined)
        continue
    end
    if cell_index < 1 || cell_index > numel(All_Cells_combined(session_index).GMM_based_classification_days)
        continue
    end

    All_Cells_combined(session_index).GMM_based_classification_days(cell_index) = GMM_based_classification(iRow);
end
end

function n_cells = infer_session_unit_count(session_entry)

candidate_fields = {'spatial_coverage_meanOFs', 'classific_firingRate', 'acg_mean', ...
    'CCGbased_CellClassfication', 'GMM_based_classification_days'};
n_cells = 0;
for iField = 1:numel(candidate_fields)
    field_name = candidate_fields{iField};
    if isfield(session_entry, field_name) && ~isempty(session_entry.(field_name))
        n_cells = max(n_cells, numel(session_entry.(field_name)));
    end
end
end

function values = get_numeric_field(session_entry, field_name, n_cells, default_value)

values = repmat(default_value, n_cells, 1);
if ~isfield(session_entry, field_name) || isempty(session_entry.(field_name))
    return
end

raw_values = double(session_entry.(field_name)(:));
n = min(n_cells, numel(raw_values));
values(1:n) = raw_values(1:n);
end

function all_cells_path = resolve_all_cells_path()

repo_root = fullfile(fileparts(mfilename('fullpath')), '..', '..');
candidate_paths = { ...
    fullfile(repo_root, 'All_Cells_combined.mat'), ...
    fullfile(repo_root, 'Data', 'All_Cells_combined.mat'), ...
    fullfile(pwd, 'All_Cells_combined.mat')};

all_cells_path = '';
for iPath = 1:numel(candidate_paths)
    if exist(candidate_paths{iPath}, 'file')
        all_cells_path = candidate_paths{iPath};
        return
    end
end

error('Could not find All_Cells_combined.mat. Pass custom_settings.allCellsPath or place it in the repository root/Data folder.')
end

function value = get_override_value(settings_struct, field_name, default_value)

value = default_value;
if nargin < 1 || isempty(settings_struct) || ~isstruct(settings_struct)
    return
end

if isfield(settings_struct, field_name) && ~isempty(settings_struct.(field_name))
    value = settings_struct.(field_name);
end
end
