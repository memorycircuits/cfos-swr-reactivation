function Classifications_Spatial_Info(custom_settings)
% Combine GMM cell type labels with place-field outputs for final classes.

if nargin < 1 || isempty(custom_settings)
    custom_settings = struct();
end

analyse_maze = char(string(get_override_value(custom_settings, 'analyseMaze', 'of')));
update_script_name = 'Classifications_Spatial_Info.m';
all_cells_path = char(string(get_override_value(custom_settings, 'allCellsPath', '')));
if isempty(all_cells_path)
    all_cells_path = resolve_all_cells_path();
end
session_info_path = char(string(get_override_value(custom_settings, 'sessionInfoPath', '')));
if isempty(session_info_path)
    session_info_path = resolve_session_info_path();
end
session_indices = get_override_value(custom_settings, 'sessionIndices', []);

loaded_all_cells = load(all_cells_path, 'All_Cells_combined');
if ~isfield(loaded_all_cells, 'All_Cells_combined')
    error('All_Cells_combined was not found in %s.', all_cells_path);
end
All_Cells_combined = loaded_all_cells.All_Cells_combined;
last_updated = load_last_updated(all_cells_path, numel(All_Cells_combined));

loaded_session_info = load(session_info_path, 'sessInfo');
if ~isfield(loaded_session_info, 'sessInfo')
    error('sessInfo was not found in %s.', session_info_path);
end
sessInfo = loaded_session_info.sessInfo;

if isempty(session_indices)
    session_indices = find(arrayfun(@(entry) isfield(entry, 'GMM_based_classification_days') && ...
        ~isempty(entry.GMM_based_classification_days), All_Cells_combined));
end

for idx = 1:numel(session_indices)
    iii = session_indices(idx);
    if iii < 1 || iii > numel(All_Cells_combined)
        warning('Skipping invalid session index: %d', iii);
        continue
    end
    if ~isfield(All_Cells_combined(iii), 'GMM_based_classification_days') || ...
            isempty(All_Cells_combined(iii).GMM_based_classification_days)
        warning('Skipping session %d because GMM_based_classification_days is missing.', iii);
        continue
    end

    placemap_dir = fullfile(sessInfo(iii).mainDir, 'processedData', 'PlaceMap');
    rate_maps = load_open_field_rate_maps(placemap_dir, analyse_maze);
    if isempty(rate_maps) || isempty(rate_maps{1})
        warning('Skipping session %d because no %s place-map file was found.', iii, analyse_maze);
        continue
    end

    GMM_based_classification_days = All_Cells_combined(iii).GMM_based_classification_days(:);
    numCells = max([numel(GMM_based_classification_days), infer_num_cells_from_rate_maps(rate_maps)]);
    GMM_based_classification_days = pad_numeric_vector(GMM_based_classification_days, numCells, 0);

    [of_place_fields, of_spatial_info, of_spatial_info_rate, of_field_size] = ...
        collect_open_field_metrics(rate_maps, numCells);

    room_ID = repmat(room_id_for_session(iii), numCells, 1);
    total_number_pfs = sum(of_place_fields, 2);
    if room_ID(1) ~= 15
        total_number_pfs = of_place_fields(:, 1);
    end

    final_cells_classification = cell(numCells, 1);
    final_classification_numeric = zeros(numCells, 1);
    place_cell_OF_combinations_code = zeros(numCells, 3);

    context_A_spatial_info = nan(numCells, 1);
    context_B_spatial_info = nan(numCells, 1);
    context_A_revisit_spatial_info = nan(numCells, 1);
    context_A_spatial_info_rate = nan(numCells, 1);
    context_B_spatial_info_rate = nan(numCells, 1);
    context_A_revisit_spatial_info_rate = nan(numCells, 1);
    context_A_field_size = nan(numCells, 1);
    context_B_field_size = nan(numCells, 1);
    context_A_revisit_field_size = nan(numCells, 1);

    for c = 1:numCells
        if GMM_based_classification_days(c) == 1
            if total_number_pfs(c) > 0
                final_cells_classification{c} = 'Place Cell';
                final_classification_numeric(c) = 2;
            else
                final_cells_classification{c} = 'Pyramidal';
                final_classification_numeric(c) = 1;
            end
        elseif GMM_based_classification_days(c) == 2
            final_cells_classification{c} = 'Interneuron';
            final_classification_numeric(c) = 3;
        else
            final_cells_classification{c} = 'Unclassified';
            final_classification_numeric(c) = 4;
        end

        if GMM_based_classification_days(c) == 1
            place_cell_OF_combinations_code(c, :) = double(of_place_fields(c, :) > 0);
        end
    end

    has_of1 = of_place_fields(:, 1) > 0;
    context_A_spatial_info(has_of1) = of_spatial_info(has_of1, 1);
    context_A_spatial_info_rate(has_of1) = of_spatial_info_rate(has_of1, 1);
    context_A_field_size(has_of1) = of_field_size(has_of1, 1);

    has_of2 = of_place_fields(:, 2) > 0;
    context_B_spatial_info(has_of2) = of_spatial_info(has_of2, 2);
    context_B_spatial_info_rate(has_of2) = of_spatial_info_rate(has_of2, 2);
    context_B_field_size(has_of2) = of_field_size(has_of2, 2);

    has_of3 = of_place_fields(:, 3) > 0;
    context_A_revisit_spatial_info(has_of3) = of_spatial_info(has_of3, 3);
    context_A_revisit_spatial_info_rate(has_of3) = of_spatial_info_rate(has_of3, 3);
    context_A_revisit_field_size(has_of3) = of_field_size(has_of3, 3);

    dataTable = table(final_cells_classification, context_A_spatial_info, context_B_spatial_info, ...
        context_A_revisit_spatial_info, context_A_spatial_info_rate, context_B_spatial_info_rate, ...
        context_A_revisit_spatial_info_rate, context_A_field_size, context_B_field_size, ...
        context_A_revisit_field_size, room_ID, ...
        'VariableNames', {'Final Classifications', 'Context A Spatial Info', ...
        'Context B Spatial Info', 'Context A Revisit Spatial Info', ...
        'Context A Spatial Info_Rate', 'Context B Spatial Info_Rate', ...
        'Context A Revisit Spatial Info_Rate', 'Context A Field Size', ...
        'Context B Field Size', 'Context A Revisit Field Size', 'room_ID'});

    putative_cells_file = fullfile(sessInfo(iii).mainDir, 'processedData', 'Classifications_Spatial_Info.mat');
    save(putative_cells_file, 'dataTable');

    All_Cells_combined(iii).final_cells_classification = final_cells_classification;
    All_Cells_combined(iii).final_classification_numeric = final_classification_numeric;
    All_Cells_combined(iii).i_number = repmat(iii, size(final_classification_numeric));
    All_Cells_combined(iii).animal = repmat(sessInfo(iii).animal, size(final_classification_numeric));
    All_Cells_combined(iii).rec_day = repmat(sessInfo(iii).day, size(final_classification_numeric));
    All_Cells_combined(iii).room_ID = room_ID;
    All_Cells_combined(iii).place_cell_OF_combinations_code = place_cell_OF_combinations_code;
    All_Cells_combined(iii).context_A_spatial_info = context_A_spatial_info;
    All_Cells_combined(iii).context_B_spatial_info = context_B_spatial_info;
    All_Cells_combined(iii).context_A_revisit_spatial_info = context_A_revisit_spatial_info;
    All_Cells_combined(iii).context_A_spatial_info_rate = context_A_spatial_info_rate;
    All_Cells_combined(iii).context_B_spatial_info_rate = context_B_spatial_info_rate;
    All_Cells_combined(iii).context_A_revisit_spatial_info_rate = context_A_revisit_spatial_info_rate;
    All_Cells_combined(iii).context_A_field_size = context_A_field_size;
    All_Cells_combined(iii).context_B_field_size = context_B_field_size;
    All_Cells_combined(iii).context_A_revisit_field_size = context_A_revisit_field_size;
    All_Cells_combined(iii).of1_place_field_numbers = of_place_fields(:, 1);
    All_Cells_combined(iii).of2_place_field_numbers = of_place_fields(:, 2);
    All_Cells_combined(iii).of3_place_field_numbers = of_place_fields(:, 3);

    last_updated = mark_classification_fields_updated(last_updated, iii, update_script_name);
end

save(all_cells_path, 'All_Cells_combined', 'last_updated', '-v7.3');
fprintf('Final cell classification fields were written to %s\n', all_cells_path);
end

function rate_maps = load_open_field_rate_maps(placemap_dir, analyse_maze)

rate_maps = cell(1, 3);
for iMap = 1:3
    map_path = fullfile(placemap_dir, sprintf('%s%d_Map.mat', analyse_maze, iMap));
    if ~exist(map_path, 'file')
        continue
    end
    loaded_map = load(map_path, 'rateMap');
    if isfield(loaded_map, 'rateMap')
        rate_maps{iMap} = loaded_map.rateMap;
    end
end
end

function [of_place_fields, of_spatial_info, of_spatial_info_rate, of_field_size] = collect_open_field_metrics(rate_maps, numCells)

of_place_fields = zeros(numCells, 3);
of_spatial_info = nan(numCells, 3);
of_spatial_info_rate = nan(numCells, 3);
of_field_size = nan(numCells, 3);

for iMap = 1:3
    rate_map = rate_maps{iMap};
    if isempty(rate_map)
        continue
    end

    of_place_fields(:, iMap) = get_rate_map_vector(rate_map, 'fieldNum', numCells, 0);
    of_spatial_info(:, iMap) = get_rate_map_vector(rate_map, 'information', numCells, NaN);
    of_spatial_info_rate(:, iMap) = get_rate_map_vector(rate_map, 'information2', numCells, NaN);
    of_field_size(:, iMap) = get_rate_map_vector(rate_map, 'sumSize', numCells, NaN);
end
end

function numCells = infer_num_cells_from_rate_maps(rate_maps)

numCells = 0;
for iMap = 1:numel(rate_maps)
    rate_map = rate_maps{iMap};
    if isempty(rate_map)
        continue
    end
    candidate_fields = {'fieldNum', 'information', 'information2', 'sumSize'};
    for iField = 1:numel(candidate_fields)
        field_name = candidate_fields{iField};
        if isfield(rate_map, field_name) && ~isempty(rate_map.(field_name))
            numCells = max(numCells, numel(rate_map.(field_name)));
        end
    end
end
end

function values = get_rate_map_vector(rate_map, field_name, numCells, default_value)

values = repmat(default_value, numCells, 1);
if ~isfield(rate_map, field_name) || isempty(rate_map.(field_name))
    return
end
raw_values = double(rate_map.(field_name)(:));
n = min(numCells, numel(raw_values));
values(1:n) = raw_values(1:n);
end

function room_id = room_id_for_session(session_index)

control_room_familiar = [143, 146, 149, 152, 163, 166, 169, 183, 186, 189, 192, 195, ...
    198, 203, 206, 209, 212, 215, 218, 223, 226, 229, 232, 235, ...
    238, 243, 246, 249, 252, 255];
control_room_novel = [144, 147, 150, 153, 164, 167, 170, 184, 187, 190, 193, 196, ...
    199, 204, 207, 210, 213, 216, 219, 224, 227, 230, 233, 236, ...
    239, 244, 247, 250, 253, 256];

if ismember(session_index, control_room_familiar)
    room_id = 19;
elseif ismember(session_index, control_room_novel)
    room_id = 20;
else
    room_id = 15;
end
end

function last_updated = load_last_updated(all_cells_path, n_sessions)

available_variables = who('-file', all_cells_path);
if any(strcmp(available_variables, 'last_updated'))
    loaded_last_updated = load(all_cells_path, 'last_updated');
    last_updated = normalize_last_updated(loaded_last_updated.last_updated, n_sessions);
else
    last_updated = repmat(struct(), n_sessions, 1);
end
end

function last_updated = normalize_last_updated(last_updated, n_sessions)

if numel(last_updated) < n_sessions
    last_updated(numel(last_updated)+1:n_sessions, 1) = struct();
end

for iSession = 1:n_sessions
    session_fields = fieldnames(last_updated(iSession));
    for iField = 1:numel(session_fields)
        field_name = session_fields{iField};
        field_value = last_updated(iSession).(field_name);

        if isstruct(field_value)
            if ~isfield(field_value, 'update_timestamp')
                field_value.update_timestamp = '';
            end
            if ~isfield(field_value, 'update_script')
                field_value.update_script = 'legacy_unknown';
            end
        else
            field_value = struct('update_timestamp', '', 'update_script', 'legacy_unknown');
        end
        last_updated(iSession).(field_name) = field_value;
    end
end
end

function last_updated = mark_classification_fields_updated(last_updated, session_index, update_script_name)

updated_fields = { ...
    'final_cells_classification', ...
    'final_classification_numeric', ...
    'i_number', ...
    'animal', ...
    'rec_day', ...
    'room_ID', ...
    'place_cell_OF_combinations_code', ...
    'context_A_spatial_info', ...
    'context_B_spatial_info', ...
    'context_A_revisit_spatial_info', ...
    'context_A_spatial_info_rate', ...
    'context_B_spatial_info_rate', ...
    'context_A_revisit_spatial_info_rate', ...
    'context_A_field_size', ...
    'context_B_field_size', ...
    'context_A_revisit_field_size', ...
    'of1_place_field_numbers', ...
    'of2_place_field_numbers', ...
    'of3_place_field_numbers' ...
    };

update_timestamp = datestr(now, 'yyyy-mm-ddTHH:MM:SS');
for iField = 1:numel(updated_fields)
    last_updated(session_index).(updated_fields{iField}) = struct( ...
        'update_timestamp', update_timestamp, ...
        'update_script', update_script_name);
end
end

function values = pad_numeric_vector(values, n, default_value)

values = double(values(:));
if numel(values) < n
    values(numel(values)+1:n, 1) = default_value;
else
    values = values(1:n);
end
end

function all_cells_path = resolve_all_cells_path()

repo_root = fullfile(fileparts(mfilename('fullpath')), '..', '..');
candidate_paths = { ...
    fullfile(repo_root, 'All_Cells_combined.mat'), ...
    fullfile(repo_root, 'Data', 'All_Cells_combined.mat'), ...
    fullfile(pwd, 'All_Cells_combined.mat')};

all_cells_path = first_existing_file(candidate_paths);
if isempty(all_cells_path)
    error('Could not find All_Cells_combined.mat. Pass custom_settings.allCellsPath or place it in the repository root/Data folder.')
end
end

function session_info_path = resolve_session_info_path()

repo_root = fullfile(fileparts(mfilename('fullpath')), '..', '..');
candidate_paths = { ...
    fullfile(repo_root, 'sessionInfo.mat'), ...
    fullfile(repo_root, 'Data', 'sessionInfo.mat'), ...
    fullfile(repo_root, 'Analysis_scripts', 'DataOrganization', 'sessionInfo.mat'), ...
    fullfile(pwd, 'sessionInfo.mat')};

session_info_path = first_existing_file(candidate_paths);
if isempty(session_info_path)
    error('Could not find sessionInfo.mat. Pass custom_settings.sessionInfoPath or place it in the repository root/Data folder.')
end
end

function path_out = first_existing_file(candidate_paths)

path_out = '';
for iPath = 1:numel(candidate_paths)
    if exist(candidate_paths{iPath}, 'file')
        path_out = candidate_paths{iPath};
        return
    end
end
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
