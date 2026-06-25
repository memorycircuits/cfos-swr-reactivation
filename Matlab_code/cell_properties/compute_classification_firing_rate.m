function result = compute_classification_firing_rate(custom_settings)
% Compute classific_firingRate for cell-type classification.

if nargin < 1 || isempty(custom_settings)
    custom_settings = struct();
end

config = load_classification_config(get_option_value(custom_settings, 'configPath', ''));
add_dependency_paths(config, custom_settings);

all_cells_path = char(string(get_option_value(custom_settings, 'allCellsPath', '')));
if isempty(all_cells_path)
    all_cells_path = resolve_all_cells_path();
end

session_info_path = char(string(get_option_value(custom_settings, 'sessionInfoPath', '')));
if isempty(session_info_path)
    session_info_path = resolve_session_info_path();
end

session_indices = get_option_value(custom_settings, 'sessionIndices', []);
sleep_folders = get_option_value(custom_settings, 'sleepFolders', {'s1', 's2'});
include_sleep = logical(get_option_value(custom_settings, 'includeSleep', true));
save_updated_all_cells = logical(get_option_value(custom_settings, 'saveUpdatedAllCells', true));
update_script_name = 'compute_classification_firing_rate.m';

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

session_summaries = repmat(empty_firing_rate_summary(), numel(session_indices), 1);

for idx = 1:numel(session_indices)
    session_index = session_indices(idx);
    session_summaries(idx).session_index = session_index;

    if session_index < 1 || session_index > numel(sessInfo) || session_index > numel(All_Cells_combined)
        session_summaries(idx).status = 'skipped_invalid_session_index';
        warning('Skipping invalid session index: %d', session_index);
        continue
    end

    main_dir = sessInfo(session_index).mainDir;
    tt_list_path = fullfile(main_dir, sessInfo(session_index).tList);
    tt_files = read_tt_list(tt_list_path);
    if isempty(tt_files)
        session_summaries(idx).status = 'skipped_empty_tetrode_list';
        warning('Session %d: no tetrode files in %s.', session_index, tt_list_path);
        continue
    end

    total_spikes = zeros(numel(tt_files), 1);
    total_duration = 0;
    processed_subsessions = {};

    open_field_folders = find_open_field_folders(main_dir);
    indata_of = load_indata_file(fullfile(main_dir, 'processedData', 'indata_of.mat'));
    for of_idx = 1:numel(open_field_folders)
        folder_name = open_field_folders{of_idx};
        [spike_times, loaded_ok] = read_subsession_spikes(fullfile(main_dir, folder_name), tt_files);
        if ~loaded_ok
            continue
        end
        total_spikes = total_spikes + count_spikes_per_cell(spike_times, numel(tt_files));
        duration = duration_from_indata(indata_of, of_idx);
        if ~isfinite(duration) || duration <= 0
            duration = infer_spike_duration(spike_times);
        end
        total_duration = total_duration + duration;
        processed_subsessions{end + 1} = folder_name; %#ok<AGROW>
    end

    if include_sleep
        indata_sleep = load_indata_file(fullfile(main_dir, 'processedData', 'indataS.mat'));
        for sleep_idx = 1:numel(sleep_folders)
            folder_name = sleep_folders{sleep_idx};
            [spike_times, loaded_ok] = read_subsession_spikes(fullfile(main_dir, folder_name), tt_files);
            if ~loaded_ok
                continue
            end
            total_spikes = total_spikes + count_spikes_per_cell(spike_times, numel(tt_files));
            duration = duration_from_indata(indata_sleep, sleep_idx);
            if ~isfinite(duration) || duration <= 0
                duration = infer_spike_duration(spike_times);
            end
            total_duration = total_duration + duration;
            processed_subsessions{end + 1} = folder_name; %#ok<AGROW>
        end
    end

    if total_duration <= 0 || ~isfinite(total_duration)
        session_summaries(idx).status = 'skipped_no_valid_duration';
        warning('Session %d: no valid recording duration found.', session_index);
        continue
    end

    firing_rate = total_spikes ./ total_duration;
    All_Cells_combined(session_index).classific_firingRate = firing_rate(:);
    last_updated = mark_fields_updated(last_updated, session_index, ...
        {'classific_firingRate'}, update_script_name);

    session_summaries(idx).status = 'processed';
    session_summaries(idx).cell_count = numel(firing_rate);
    session_summaries(idx).total_duration_seconds = total_duration;
    session_summaries(idx).processed_subsessions = processed_subsessions;
end

if save_updated_all_cells
    save_updated_all_cells_file(all_cells_path, All_Cells_combined, last_updated);
end

result = struct( ...
    'allCellsPath', all_cells_path, ...
    'sessionInfoPath', session_info_path, ...
    'includeSleep', include_sleep, ...
    'sessions', session_summaries);

fprintf('Computed classific_firingRate for %d sessions.\n', ...
    nnz(strcmp({session_summaries.status}, 'processed')));
end

function [spike_times, loaded_ok] = read_subsession_spikes(session_dir, tt_files)

spike_times = {};
loaded_ok = false;
if exist(session_dir, 'dir') ~= 7
    warning('Skipping missing subsession folder: %s', session_dir);
    return
end

try
    spike_data = readSpikeDataOnly(session_dir, tt_files);
    spike_times = fixSpikes(spike_data);
catch ME
    warning('Could not read spikes from %s: %s', session_dir, ME.message);
    return
end

spike_times = spike_times(:);
loaded_ok = true;
end

function counts = count_spikes_per_cell(spike_times, cell_count)

counts = zeros(cell_count, 1);
for cell_idx = 1:min(cell_count, numel(spike_times))
    counts(cell_idx) = numel(spike_times{cell_idx});
end
end

function duration = duration_from_indata(indata, index)

duration = NaN;
if isempty(indata) || index > numel(indata) || ~isfield(indata(index), 't') || isempty(indata(index).t)
    return
end

t = double(indata(index).t(:));
t = t(isfinite(t));
if numel(t) >= 2
    duration = max(t) - min(t);
end
end

function duration = infer_spike_duration(spike_times)

all_times = [];
for cell_idx = 1:numel(spike_times)
    all_times = [all_times; spike_times{cell_idx}(:)]; %#ok<AGROW>
end
all_times = double(all_times(isfinite(all_times)));
if numel(all_times) >= 2
    duration = max(all_times) - min(all_times);
else
    duration = NaN;
end
end

function indata = load_indata_file(indata_path)

indata = [];
if exist(indata_path, 'file') ~= 2
    return
end

loaded = load(indata_path, 'indata');
if isfield(loaded, 'indata')
    indata = loaded.indata;
end
end

function folder_names = find_open_field_folders(main_dir)

folder_names = {};
listing = dir(main_dir);
for entry_idx = 1:numel(listing)
    if ~listing(entry_idx).isdir
        continue
    end
    tokens = regexp(listing(entry_idx).name, '^of(\d+)$', 'tokens', 'once');
    if isempty(tokens)
        continue
    end
    folder_names{end + 1, 1} = listing(entry_idx).name; %#ok<AGROW>
end

folder_numbers = cellfun(@(name) str2double(regexp(name, '\d+', 'match', 'once')), folder_names);
[~, order] = sort(folder_numbers);
folder_names = folder_names(order);
end

function tt_files = read_tt_list(tt_list_path)

tt_files = {};
fid = fopen(tt_list_path);
if fid < 0
    warning('Could not open tetrode list: %s', tt_list_path);
    return
end

cleanup = onCleanup(@() fclose(fid));
while true
    line_value = fgetl(fid);
    if ~ischar(line_value)
        break
    end
    line_value = strtrim(line_value);
    if isempty(line_value)
        continue
    end
    tt_files{end + 1, 1} = line_value; %#ok<AGROW>
end
delete(cleanup);
end

function summary = empty_firing_rate_summary()

summary = struct( ...
    'session_index', NaN, ...
    'status', '', ...
    'cell_count', 0, ...
    'total_duration_seconds', NaN, ...
    'processed_subsessions', {{}} ...
    );
end

function config = load_classification_config(config_path)

if nargin < 1 || isempty(config_path)
    config_path = fullfile(fileparts(mfilename('fullpath')), '..', 'classification', 'classification_config.json');
end

config = struct('cellExplorerPath', '', 'mclustPath', '');
if exist(config_path, 'file') ~= 2
    return
end

try
    decoded = jsondecode(fileread(config_path));
catch ME
    warning('Could not read classification config %s: %s', config_path, ME.message);
    return
end

if isfield(decoded, 'cellExplorerPath')
    config.cellExplorerPath = char(string(decoded.cellExplorerPath));
end
if isfield(decoded, 'mclustPath')
    config.mclustPath = char(string(decoded.mclustPath));
end
end

function add_dependency_paths(config, custom_settings)

matlab_code_folder = fullfile(fileparts(mfilename('fullpath')), '..');
if exist(matlab_code_folder, 'dir') == 7
    addpath(genpath(matlab_code_folder));
end

mclust_path = char(string(get_option_value(custom_settings, 'mclustPath', config.mclustPath)));
add_dependency_path(mclust_path, 'MClust');

additional_paths = get_option_value(custom_settings, 'additionalPaths', {});
if ischar(additional_paths) || isstring(additional_paths)
    additional_paths = cellstr(string(additional_paths));
end
for path_idx = 1:numel(additional_paths)
    add_dependency_path(additional_paths{path_idx}, 'additional dependency');
end
end

function add_dependency_path(path_value, dependency_name)

path_value = char(string(path_value));
if isempty(path_value)
    return
end

if exist(path_value, 'dir') == 7
    addpath(genpath(path_value));
else
    warning('%s path does not exist: %s', dependency_name, path_value);
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
