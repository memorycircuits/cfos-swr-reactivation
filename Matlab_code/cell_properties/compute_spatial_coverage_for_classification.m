function result = compute_spatial_coverage_for_classification(custom_settings)
% Compute spatial coverage fields used by the GMM classifier.

if nargin < 1 || isempty(custom_settings)
    custom_settings = struct();
end

add_cell_properties_path();

all_cells_path = char(string(get_option_value(custom_settings, 'allCellsPath', '')));
if isempty(all_cells_path)
    all_cells_path = resolve_all_cells_path();
end

session_info_path = char(string(get_option_value(custom_settings, 'sessionInfoPath', '')));
if isempty(session_info_path)
    session_info_path = resolve_session_info_path();
end

maze_prefix = char(string(get_option_value(custom_settings, 'mazePrefix', 'of')));
session_indices = get_option_value(custom_settings, 'sessionIndices', []);
threshold = double(get_option_value(custom_settings, 'threshold', 0.75));
save_updated_all_cells = logical(get_option_value(custom_settings, 'saveUpdatedAllCells', true));
update_script_name = 'compute_spatial_coverage_for_classification.m';

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
    session_indices = 1:numel(sessInfo);
end

session_summaries = repmat(empty_spatial_summary(), numel(session_indices), 1);

for idx = 1:numel(session_indices)
    session_index = session_indices(idx);
    session_summaries(idx).session_index = session_index;

    if session_index < 1 || session_index > numel(sessInfo) || session_index > numel(All_Cells_combined)
        session_summaries(idx).status = 'skipped_invalid_session_index';
        warning('Skipping invalid session index: %d', session_index);
        continue
    end

    main_dir = sessInfo(session_index).mainDir;
    place_map_dir = fullfile(main_dir, 'processedData', 'PlaceMap');
    map_files = find_open_field_map_files(place_map_dir, maze_prefix);
    if isempty(map_files)
        session_summaries(idx).status = 'skipped_no_place_maps';
        warning('Session %d: no %s*_Map.mat files found in %s.', ...
            session_index, maze_prefix, place_map_dir);
        continue
    end

    [All_Cells_combined, written_fields, coverage_matrix] = ...
        write_session_spatial_coverage(All_Cells_combined, session_index, map_files, threshold);

    All_Cells_combined(session_index).spatial_coverage_meanOFs = ...
        row_nanmean(coverage_matrix);
    written_fields{end + 1} = 'spatial_coverage_meanOFs';

    last_updated = mark_fields_updated(last_updated, session_index, ...
        written_fields, update_script_name);

    session_summaries(idx).status = 'processed';
    session_summaries(idx).map_count = numel(map_files);
    session_summaries(idx).cell_count = numel(All_Cells_combined(session_index).spatial_coverage_meanOFs);
    session_summaries(idx).written_fields = written_fields;
end

if save_updated_all_cells
    save_updated_all_cells_file(all_cells_path, All_Cells_combined, last_updated);
end

result = struct( ...
    'allCellsPath', all_cells_path, ...
    'sessionInfoPath', session_info_path, ...
    'mazePrefix', maze_prefix, ...
    'threshold', threshold, ...
    'sessions', session_summaries);

fprintf('Computed spatial_coverage_meanOFs for %d sessions.\n', ...
    nnz(strcmp({session_summaries.status}, 'processed')));
end

function [All_Cells_combined, written_fields, coverage_matrix] = ...
        write_session_spatial_coverage(All_Cells_combined, session_index, map_files, threshold)

coverage_columns = {};
written_fields = {};
cell_count = 0;

for map_idx = 1:numel(map_files)
    loaded = load(map_files(map_idx).path);
    if ~isfield(loaded, 'rateMap')
        warning('Skipping map file without rateMap variable: %s', map_files(map_idx).path);
        continue
    end

    rate_map = loaded.rateMap;
    maps = extract_rate_maps(rate_map);
    coverage_values = NaN(numel(maps), 1);
    for cell_idx = 1:numel(maps)
        coverage_values(cell_idx) = computeSpatialCoverage(maps{cell_idx}, ...
            'Threshold', threshold, ...
            'DoPlot', false);
    end

    cell_count = max(cell_count, numel(coverage_values));
    field_name = sprintf('spatial_coverage_%s%d', map_files(map_idx).prefix, map_files(map_idx).number);
    All_Cells_combined(session_index).(field_name) = coverage_values(:);
    written_fields{end + 1} = field_name; %#ok<AGROW>
    coverage_columns{end + 1} = coverage_values(:); %#ok<AGROW>

    if isfield(rate_map, 'tFileNames') && ~isempty(rate_map.tFileNames)
        All_Cells_combined(session_index).tFileNames = rate_map.tFileNames(:);
        if ~any(strcmp(written_fields, 'tFileNames'))
            written_fields{end + 1} = 'tFileNames'; %#ok<AGROW>
        end
    end
end

coverage_matrix = NaN(cell_count, numel(coverage_columns));
for col_idx = 1:numel(coverage_columns)
    values = coverage_columns{col_idx};
    coverage_matrix(1:numel(values), col_idx) = values;
end
end

function maps = extract_rate_maps(rate_map)

if ~isstruct(rate_map) || ~isfield(rate_map, 'map')
    error('rateMap must be a struct with a map field.');
end

map_input = rate_map.map;
if iscell(map_input)
    maps = map_input(:);
elseif isnumeric(map_input)
    if ismatrix(map_input)
        maps = {map_input};
    elseif ndims(map_input) == 3
        maps = cell(size(map_input, 3), 1);
        for map_idx = 1:size(map_input, 3)
            maps{map_idx} = map_input(:, :, map_idx);
        end
    else
        error('Numeric rate maps must be 2-D or 3-D.');
    end
else
    error('Unsupported rateMap.map type.');
end
end

function map_files = find_open_field_map_files(place_map_dir, maze_prefix)

map_files = struct('path', {}, 'prefix', {}, 'number', {});
if exist(place_map_dir, 'dir') ~= 7
    return
end

listing = dir(fullfile(place_map_dir, sprintf('%s*_Map.mat', maze_prefix)));
for file_idx = 1:numel(listing)
    tokens = regexp(listing(file_idx).name, sprintf('^(%s)(\\d+)_Map\\.mat$', regexptranslate('escape', maze_prefix)), 'tokens', 'once');
    if isempty(tokens)
        continue
    end
    map_files(end + 1).path = fullfile(listing(file_idx).folder, listing(file_idx).name); %#ok<AGROW>
    map_files(end).prefix = tokens{1};
    map_files(end).number = str2double(tokens{2});
end

[~, order] = sort([map_files.number]);
map_files = map_files(order);
end

function values = row_nanmean(matrix_values)

if isempty(matrix_values)
    values = [];
    return
end

valid_counts = sum(isfinite(matrix_values), 2);
matrix_values(~isfinite(matrix_values)) = 0;
values = sum(matrix_values, 2) ./ valid_counts;
values(valid_counts == 0) = NaN;
end

function summary = empty_spatial_summary()

summary = struct( ...
    'session_index', NaN, ...
    'status', '', ...
    'map_count', 0, ...
    'cell_count', 0, ...
    'written_fields', {{}} ...
    );
end

function add_cell_properties_path()

cell_properties_folder = fileparts(mfilename('fullpath'));
if exist(cell_properties_folder, 'dir')
    addpath(cell_properties_folder);
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
    error('Could not find All_Cells_combined.mat. Pass custom_settings.allCellsPath or place it in the repository root/Data folder.');
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
    error('Could not find sessionInfo.mat. Pass custom_settings.sessionInfoPath or place it in the repository root/Data folder.');
end
end

function path_out = first_existing_file(candidate_paths)

path_out = '';
for path_idx = 1:numel(candidate_paths)
    if exist(candidate_paths{path_idx}, 'file') == 2
        path_out = candidate_paths{path_idx};
        return
    end
end
end

function last_updated = load_last_updated(all_cells_path, session_count)

last_updated = struct();
if exist(all_cells_path, 'file') ~= 2
    return
end

available_variables = who('-file', all_cells_path);
if any(strcmp(available_variables, 'last_updated'))
    loaded = load(all_cells_path, 'last_updated');
    if isstruct(loaded.last_updated)
        last_updated = loaded.last_updated;
    end
end

last_updated = ensure_last_updated_session_capacity(last_updated, session_count);
end

function save_updated_all_cells_file(all_cells_path, All_Cells_combined, last_updated)

if exist(all_cells_path, 'file') == 2
    save(all_cells_path, 'All_Cells_combined', 'last_updated', '-append');
else
    save(all_cells_path, 'All_Cells_combined', 'last_updated', '-v7.3');
end
end

function last_updated = mark_fields_updated(last_updated, session_index, fields, update_script_name)

last_updated = ensure_last_updated_session_capacity(last_updated, session_index);
timestamp = datestr(now, 'yyyy-mm-ddTHH:MM:SS');
for field_idx = 1:numel(fields)
    field_name = fields{field_idx};
    if ~isfield(last_updated.sessions(session_index).fields, field_name)
        last_updated.sessions(session_index).fields.(field_name) = struct();
    end
    last_updated.sessions(session_index).fields.(field_name).timestamp = timestamp;
    last_updated.sessions(session_index).fields.(field_name).script = update_script_name;
end
end

function last_updated = ensure_last_updated_session_capacity(last_updated, session_count)

if ~isfield(last_updated, 'sessions') || numel(last_updated.sessions) < session_count
    if ~isfield(last_updated, 'sessions')
        last_updated.sessions = struct('fields', {});
    end
    current_count = numel(last_updated.sessions);
    for session_idx = current_count + 1:session_count
        last_updated.sessions(session_idx).fields = struct();
    end
end
end

function value = get_option_value(settings_struct, field_name, default_value)

value = default_value;
if nargin < 1 || isempty(settings_struct) || ~isstruct(settings_struct)
    return
end

if isfield(settings_struct, field_name) && ~isempty(settings_struct.(field_name))
    value = settings_struct.(field_name);
end
end
