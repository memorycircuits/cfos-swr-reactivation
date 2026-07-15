function manifest = collect_swr_cell_metrics_demo_data(session_info_path, session_index, target_root, varargin)
%COLLECT_SWR_CELL_METRICS_DEMO_DATA Package one sleep session for the SWR metrics demo.
%   MANIFEST = COLLECT_SWR_CELL_METRICS_DEMO_DATA(SESSION_INFO_PATH,
%   SESSION_INDEX, TARGET_ROOT) copies the minimal S1 inputs required by
%   compute_swr_cell_metrics into TARGET_ROOT/swr_cell_metrics_demo.
%
%   Name-value options:
%     SleepFolder       Sleep folder to package. Default: 's1'.
%     UnitIndices       Indices from the source tList to include. Default: all.
%     MaxCSCSizeMB      Maximum permitted CSC-file size. Default: 100.
%     Overwrite         Replace an existing output folder. Default: false.

parser = inputParser;
addParameter(parser, 'SleepFolder', 's1', @(x) ischar(x) || isstring(x));
addParameter(parser, 'UnitIndices', [], @(x) isnumeric(x) && isvector(x) && all(x >= 1));
addParameter(parser, 'MaxCSCSizeMB', 100, @(x) isnumeric(x) && isscalar(x) && x > 0);
addParameter(parser, 'Overwrite', false, @(x) islogical(x) && isscalar(x));
parse(parser, varargin{:});
options = parser.Results;

session_info_path = char(string(session_info_path));
target_root = char(string(target_root));
sleep_folder = char(string(options.SleepFolder));
unit_indices = unique(double(options.UnitIndices(:)'));

if exist(session_info_path, 'file') ~= 2
    error('SWRDemoCollector:MissingSessionInfo', 'sessionInfo file not found: %s', session_info_path);
end

loaded_session_info = load(session_info_path, 'sessInfo');
if ~isfield(loaded_session_info, 'sessInfo')
    error('SWRDemoCollector:MissingSessInfo', 'sessInfo was not found in %s.', session_info_path);
end
sessInfo = loaded_session_info.sessInfo;

if session_index < 1 || session_index > numel(sessInfo) || session_index ~= floor(session_index)
    error('SWRDemoCollector:InvalidSessionIndex', ...
        'session_index must be an integer between 1 and %d.', numel(sessInfo));
end

source_session = sessInfo(session_index);
required_fields = {'mainDir', 'tList', 'cellLayerChann'};
for field_index = 1:numel(required_fields)
    field_name = required_fields{field_index};
    if ~isfield(source_session, field_name) || isempty(source_session.(field_name))
        error('SWRDemoCollector:MissingSessionField', ...
            'sessInfo(%d).%s is required.', session_index, field_name);
    end
end

main_dir = char(string(source_session.mainDir));
t_list_name = char(string(source_session.tList));
csc_channel = double(source_session.cellLayerChann);
if exist(main_dir, 'dir') ~= 7
    error('SWRDemoCollector:MissingMainDirectory', 'Session directory not found: %s', main_dir);
end

source_sleep_dir = fullfile(main_dir, sleep_folder);
if exist(source_sleep_dir, 'dir') ~= 7
    error('SWRDemoCollector:MissingSleepDirectory', 'Sleep directory not found: %s', source_sleep_dir);
end

source_t_list = fullfile(main_dir, t_list_name);
tt_files = read_tt_list(source_t_list);
if isempty(tt_files)
    error('SWRDemoCollector:EmptyTList', 'No units found in %s.', source_t_list);
end
if isempty(unit_indices)
    unit_indices = 1:numel(tt_files);
end
if any(unit_indices > numel(tt_files))
    error('SWRDemoCollector:InvalidUnitIndex', ...
        'UnitIndices must be between 1 and %d for %s.', numel(tt_files), source_t_list);
end

selected_tt_files = tt_files(unit_indices);
source_spike_files = cell(numel(selected_tt_files), 1);
destination_spike_names = cell(numel(selected_tt_files), 1);
for unit_index = 1:numel(selected_tt_files)
    source_spike_files{unit_index} = resolve_spike_file(source_sleep_dir, main_dir, selected_tt_files{unit_index});
    [~, file_name, extension] = fileparts(source_spike_files{unit_index});
    destination_spike_names{unit_index} = [file_name extension];
end
if numel(unique(destination_spike_names)) ~= numel(destination_spike_names)
    error('SWRDemoCollector:DuplicateSpikeNames', ...
        'Selected spike files have duplicate filenames after packaging. Choose different UnitIndices.');
end

source_csc_file = fullfile(source_sleep_dir, sprintf('CSC%d.ncs', csc_channel));
if exist(source_csc_file, 'file') ~= 2
    error('SWRDemoCollector:MissingCSC', 'CSC file not found: %s', source_csc_file);
end
csc_info = dir(source_csc_file);
if csc_info.bytes > options.MaxCSCSizeMB * 1024^2
    error('SWRDemoCollector:CSCTooLarge', ...
        ['CSC file is %.1f MB, above the MaxCSCSizeMB limit of %.1f MB. ', ...
        'Choose a shorter source recording or raise the limit explicitly.'], ...
        csc_info.bytes / 1024^2, options.MaxCSCSizeMB);
end

source_indata_file = fullfile(main_dir, 'processedData', 'indataS.mat');
indata = load_sleep_time_vector(source_indata_file, source_session, sleep_folder);

source_swr_file = resolve_swr_file(source_sleep_dir);
slow_cHFOs = load_and_sanitize_slow_events(source_swr_file);
if isempty(slow_cHFOs)
    error('SWRDemoCollector:NoSlowEvents', 'No valid slow_cHFOs were found in %s.', source_swr_file);
end

demo_root = fullfile(target_root, 'swr_cell_metrics_demo');
if exist(demo_root, 'dir') == 7
    if ~options.Overwrite
        error('SWRDemoCollector:OutputExists', ...
            'Output folder already exists: %s. Set Overwrite to true to replace it.', demo_root);
    end
    rmdir(demo_root, 's');
end

demo_session_dir = fullfile(demo_root, 'session');
demo_sleep_dir = fullfile(demo_session_dir, sleep_folder);
demo_processed_dir = fullfile(demo_session_dir, 'processedData');
demo_sleep_processed_dir = fullfile(demo_sleep_dir, 'processedData');
mkdir(demo_sleep_processed_dir);
mkdir(demo_processed_dir);

copyfile(source_csc_file, fullfile(demo_sleep_dir, sprintf('CSC%d.ncs', csc_channel)));
for unit_index = 1:numel(source_spike_files)
    copyfile(source_spike_files{unit_index}, fullfile(demo_sleep_dir, destination_spike_names{unit_index}));
end
write_tt_list(fullfile(demo_session_dir, 'tList.txt'), destination_spike_names);

save(fullfile(demo_processed_dir, 'indataS.mat'), 'indata', '-v7');
save(fullfile(demo_sleep_processed_dir, '_allE_numSD3.5_HighPwrCycles4.mat'), 'slow_cHFOs', '-v7');

sessInfo = struct( ...
    'mainDir', '__DEMO_SESSION_ROOT__', ...
    'tList', 'tList.txt', ...
    'cellLayerChann', csc_channel);
save(fullfile(demo_root, 'sessionInfo_template.mat'), 'sessInfo', '-v7');

All_Cells_combined = struct();
last_updated = struct();
save(fullfile(demo_root, 'All_Cells_combined.mat'), 'All_Cells_combined', 'last_updated', '-v7');

manifest = struct();
manifest.source_session_index = session_index;
manifest.sleep_folder = sleep_folder;
manifest.csc_channel = csc_channel;
manifest.selected_unit_indices = unit_indices;
manifest.selected_unit_files = destination_spike_names;
manifest.slow_event_count = numel(slow_cHFOs);
manifest.first_slow_event_start = slow_cHFOs(1).start_ts;
manifest.last_slow_event_stop = slow_cHFOs(end).stop_ts;
manifest.indata_start = indata(1).t(1);
manifest.indata_stop = indata(1).t(end);
manifest.csc_size_bytes = csc_info.bytes;
manifest.output_folder = 'swr_cell_metrics_demo';
manifest.required_runtime_dependencies = {'MATLAB Signal Processing Toolbox'};
save(fullfile(demo_root, 'manifest.mat'), 'manifest', '-v7');
write_manifest_text(fullfile(demo_root, 'manifest.txt'), manifest);

fprintf('Packaged SWR cell-metrics demo data in:\n%s\n', demo_root);
fprintf('Included %d slow events and %d units.\n', ...
    manifest.slow_event_count, numel(manifest.selected_unit_files));
end

function tt_files = read_tt_list(t_list_path)

if exist(t_list_path, 'file') ~= 2
    error('SWRDemoCollector:MissingTList', 'tList file not found: %s', t_list_path);
end

fid = fopen(t_list_path, 'r');
if fid < 0
    error('SWRDemoCollector:UnreadableTList', 'Could not open tList file: %s', t_list_path);
end
cleanup = onCleanup(@() fclose(fid));

tt_files = {};
while true
    line_value = fgetl(fid);
    if ~ischar(line_value)
        break
    end
    line_value = strtrim(line_value);
    if ~isempty(line_value)
        tt_files{end + 1, 1} = line_value; %#ok<AGROW>
    end
end
end

function spike_file = resolve_spike_file(sleep_dir, main_dir, t_list_entry)

t_list_entry = char(string(t_list_entry));
candidates = {t_list_entry, fullfile(sleep_dir, t_list_entry), fullfile(main_dir, t_list_entry)};
for candidate_index = 1:numel(candidates)
    candidate = candidates{candidate_index};
    if exist(candidate, 'file') == 2
        spike_file = candidate;
        return
    end
end

error('SWRDemoCollector:MissingSpikeFile', ...
    'Could not locate the spike file listed in tList: %s', t_list_entry);
end

function indata = load_sleep_time_vector(indata_file, source_session, sleep_folder)

if exist(indata_file, 'file') ~= 2
    error('SWRDemoCollector:MissingInData', 'indataS file not found: %s', indata_file);
end

loaded = load(indata_file, 'indata');
if ~isfield(loaded, 'indata') || isempty(loaded.indata)
    error('SWRDemoCollector:MissingInDataVariable', 'indata was not found in %s.', indata_file);
end

sleep_index = resolve_sleep_index(source_session, sleep_folder, numel(loaded.indata));
if sleep_index > numel(loaded.indata) || ~isfield(loaded.indata(sleep_index), 't') || ...
        isempty(loaded.indata(sleep_index).t)
    error('SWRDemoCollector:MissingSleepTime', ...
        'indata(%d).t is required for the selected sleep folder.', sleep_index);
end

time_vector = double(loaded.indata(sleep_index).t(:));
if numel(time_vector) < 2 || any(~isfinite(time_vector)) || time_vector(end) <= time_vector(1)
    error('SWRDemoCollector:InvalidSleepTime', 'indata(%d).t must be a finite increasing time vector.', sleep_index);
end

indata = struct('t', time_vector);
end

function sleep_index = resolve_sleep_index(source_session, sleep_folder, indata_count)

sleep_index = NaN;
if isfield(source_session, 'sleepDirs') && ~isempty(source_session.sleepDirs)
    sleep_dirs = cellstr(string(source_session.sleepDirs));
    found_index = find(strcmpi(sleep_dirs, sleep_folder), 1);
    if ~isempty(found_index)
        sleep_index = found_index;
    end
end

if isnan(sleep_index)
    tokens = regexp(lower(sleep_folder), '^s(\d+)$', 'tokens', 'once');
    if ~isempty(tokens)
        sleep_index = str2double(tokens{1});
    end
end

if isnan(sleep_index) || sleep_index < 1 || sleep_index > indata_count
    error('SWRDemoCollector:UnknownSleepFolder', ...
        'Could not resolve %s to an indata entry.', sleep_folder);
end
end

function swr_file = resolve_swr_file(sleep_dir)

processed_dir = fullfile(sleep_dir, 'processedData');
expected_file = fullfile(processed_dir, '_allE_numSD3.5_HighPwrCycles4.mat');
if is_valid_swr_file(expected_file)
    swr_file = expected_file;
    return
end

listing = dir(fullfile(processed_dir, '_allE_numSD*_HighPwrCycles*.mat'));
valid_files = {};
for file_index = 1:numel(listing)
    candidate = fullfile(listing(file_index).folder, listing(file_index).name);
    if is_valid_swr_file(candidate)
        valid_files{end + 1} = candidate; %#ok<AGROW>
    end
end

if isempty(valid_files)
    error('SWRDemoCollector:MissingSWRFile', ...
        'No SWR event file containing slow_cHFOs was found in %s.', processed_dir);
end
if numel(valid_files) > 1
    error('SWRDemoCollector:AmbiguousSWRFile', ...
        'Multiple SWR event files containing slow_cHFOs were found in %s.', processed_dir);
end
swr_file = valid_files{1};
end

function tf = is_valid_swr_file(file_path)

tf = exist(file_path, 'file') == 2;
if ~tf
    return
end

variables = who('-file', file_path);
tf = any(strcmp(variables, 'slow_cHFOs'));
end

function slow_cHFOs = load_and_sanitize_slow_events(swr_file)

loaded = load(swr_file, 'slow_cHFOs');
if ~isfield(loaded, 'slow_cHFOs') || ~isstruct(loaded.slow_cHFOs)
    error('SWRDemoCollector:InvalidSlowEvents', 'slow_cHFOs is missing or invalid in %s.', swr_file);
end

source_events = loaded.slow_cHFOs(:)';
slow_cHFOs = repmat(struct('start_ts', NaN, 'stop_ts', NaN, 'duration', NaN, 'peaks_sd', NaN), ...
    1, numel(source_events));
for event_index = 1:numel(source_events)
    event = source_events(event_index);
    start_ts = get_scalar_event_value(event, 'start_ts', NaN);
    stop_ts = get_scalar_event_value(event, 'stop_ts', NaN);
    if ~isfinite(start_ts) || ~isfinite(stop_ts) || stop_ts <= start_ts
        error('SWRDemoCollector:InvalidSlowEvent', ...
            'slow_cHFOs(%d) does not have valid start_ts and stop_ts values.', event_index);
    end

    duration = get_scalar_event_value(event, 'duration', stop_ts - start_ts);
    if ~isfinite(duration) || duration <= 0
        duration = stop_ts - start_ts;
    end

    slow_cHFOs(event_index).start_ts = start_ts;
    slow_cHFOs(event_index).stop_ts = stop_ts;
    slow_cHFOs(event_index).duration = duration;
    slow_cHFOs(event_index).peaks_sd = get_scalar_event_value(event, 'peaks_sd', NaN);
end

[~, sort_index] = sort([slow_cHFOs.start_ts]);
slow_cHFOs = slow_cHFOs(sort_index);
end

function value = get_scalar_event_value(event, field_name, default_value)

value = default_value;
if isfield(event, field_name) && ~isempty(event.(field_name))
    candidate = double(event.(field_name));
    if isscalar(candidate)
        value = candidate;
    end
end
end

function write_tt_list(output_file, tt_files)

fid = fopen(output_file, 'w');
if fid < 0
    error('SWRDemoCollector:WriteTListFailed', 'Could not create %s.', output_file);
end
cleanup = onCleanup(@() fclose(fid));
for file_index = 1:numel(tt_files)
    fprintf(fid, '%s\n', tt_files{file_index});
end
end

function write_manifest_text(output_file, manifest)

fid = fopen(output_file, 'w');
if fid < 0
    error('SWRDemoCollector:WriteManifestFailed', 'Could not create %s.', output_file);
end
cleanup = onCleanup(@() fclose(fid));

fprintf(fid, 'SWR cell-metrics demo manifest\n');
fprintf(fid, 'Source session index: %d\n', manifest.source_session_index);
fprintf(fid, 'Sleep folder: %s\n', manifest.sleep_folder);
fprintf(fid, 'CSC channel: %d\n', manifest.csc_channel);
fprintf(fid, 'Selected unit indices: %s\n', mat2str(manifest.selected_unit_indices));
fprintf(fid, 'Selected unit files: %s\n', strjoin(manifest.selected_unit_files, ', '));
fprintf(fid, 'Slow event count: %d\n', manifest.slow_event_count);
fprintf(fid, 'First slow-event start: %.12g\n', manifest.first_slow_event_start);
fprintf(fid, 'Last slow-event stop: %.12g\n', manifest.last_slow_event_stop);
fprintf(fid, 'CSC size: %.2f MB\n', manifest.csc_size_bytes / 1024^2);
fprintf(fid, 'Output folder: %s\n', manifest.output_folder);
end
