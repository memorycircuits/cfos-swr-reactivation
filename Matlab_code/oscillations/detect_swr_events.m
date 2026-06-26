function result = detect_swr_events(custom_settings)
% Detect sleep SWR/HFO events and write per-session event files.

if nargin < 1 || isempty(custom_settings)
    custom_settings = struct();
end

config = load_oscillation_config(get_option_value(custom_settings, 'configPath', ''));
add_dependency_paths(config, custom_settings);

session_info_path = char(string(get_option_value(custom_settings, 'sessionInfoPath', '')));
if isempty(session_info_path)
    session_info_path = resolve_session_info_path();
end

loaded_session_info = load(session_info_path, 'sessInfo');
if ~isfield(loaded_session_info, 'sessInfo')
    error('sessInfo was not found in %s.', session_info_path);
end
sessInfo = loaded_session_info.sessInfo;

session_indices = get_option_value(custom_settings, 'sessionIndices', []);
if isempty(session_indices)
    session_indices = 1:numel(sessInfo);
end

sleep_folders = get_option_value(custom_settings, 'sleepFolders', {'s1', 's2'});
sleep_folders = normalize_text_cell(sleep_folders);
detection_threshold = resolve_detection_threshold(custom_settings);
overwrite_existing = logical(get_option_value(custom_settings, 'overwriteExisting', false));
save_hfo_rate_file = logical(get_option_value(custom_settings, 'saveHFORateFile', true));

session_summaries = repmat(empty_detection_summary(), numel(session_indices), 1);

for session_counter = 1:numel(session_indices)
    session_index = session_indices(session_counter);
    session_summaries(session_counter).session_index = session_index;

    if session_index < 1 || session_index > numel(sessInfo)
        session_summaries(session_counter).status = 'skipped_invalid_session_index';
        warning('Skipping invalid session index: %d', session_index);
        continue
    end

    sleep_summaries = repmat(empty_sleep_detection_summary(), numel(sleep_folders), 1);
    for sleep_idx = 1:numel(sleep_folders)
        sleep_label = char(string(sleep_folders{sleep_idx}));
        sleep_summaries(sleep_idx) = detect_session_sleep_swr_events( ...
            sessInfo(session_index), sleep_label, sleep_idx, detection_threshold, ...
            overwrite_existing, save_hfo_rate_file);
    end

    session_summaries(session_counter).status = summarize_session_status(sleep_summaries);
    session_summaries(session_counter).sleep = sleep_summaries;
end

result = struct( ...
    'sessionInfoPath', session_info_path, ...
    'detectionThreshold', detection_threshold, ...
    'sessions', session_summaries);

fprintf('Detected SWR events for %d sessions.\n', nnz(strcmp({session_summaries.status}, 'processed')));
end

function summary = detect_session_sleep_swr_events(sess_entry, sleep_label, sleep_idx, detection_threshold, overwrite_existing, save_hfo_rate_file)

summary = empty_sleep_detection_summary();
summary.sleep_label = sleep_label;

processed_folder = fullfile(sess_entry.mainDir, sleep_label, 'processedData');
if exist(processed_folder, 'dir') ~= 7
    mkdir(processed_folder);
end

swr_file_name = build_swr_event_filename(detection_threshold);
swr_file = fullfile(processed_folder, swr_file_name);
summary.swr_file = swr_file;

if exist(swr_file, 'file') == 2 && ~overwrite_existing
    summary.status = 'skipped_existing_swr_file';
    loaded = load(swr_file, 'slow_cHFOs');
    if isfield(loaded, 'slow_cHFOs')
        summary.slow_event_count = numel(loaded.slow_cHFOs);
    end
    return
end

indata_file = fullfile(sess_entry.mainDir, 'processedData', 'indataS.mat');
if exist(indata_file, 'file') ~= 2
    summary.status = 'skipped_missing_indataS';
    warning('Missing indataS file: %s', indata_file);
    return
end

try
    loaded_indata = load(indata_file, 'indata');
    indata_idx = resolve_sleep_indata_index(sleep_label, sleep_idx);
    if ~isfield(loaded_indata, 'indata') || isnan(indata_idx) || indata_idx < 1 || indata_idx > numel(loaded_indata.indata)
        summary.status = 'skipped_invalid_indataS';
        warning('Could not load %s indata index %d from %s.', sleep_label, indata_idx, indata_file);
        return
    end
    indata = loaded_indata.indata(indata_idx);

    lfp = load_sleep_lfp(sess_entry.mainDir, sleep_label, sess_entry.cellLayerChann);
    [lfp.ts, lfp.samp] = align_lfp_with_indata(indata.t, lfp.ts, lfp.samp);

    velocities = resize_to_lfp_samples(indata.v, numel(lfp.samp));
    quiet_mask = velocities < 2;
    fast_mask = velocities > 2;

    Hd = detectHFOs.blanco_bp(lfp.sampFreq);
    filter_numerator = Hd.Numerator;
    [slow_cHFOs, fast_cHFOs, all_cHFOs, lfp_bp] = detect_swr_events_blanco( ...
        lfp, filter_numerator, detection_threshold, velocities);
catch ME
    summary.status = 'failed';
    summary.error_message = ME.message;
    warning('SWR detection failed for %s %s: %s', sess_entry.mainDir, sleep_label, ME.message);
    return
end

time_spent = numel(lfp.samp) / lfp.sampFreq;
quiet_time_spent = sum(quiet_mask) / lfp.sampFreq;
fast_time_spent = sum(fast_mask) / lfp.sampFreq;
quiet_percent = quiet_time_spent / time_spent;
fast_percent = fast_time_spent / time_spent;
if quiet_time_spent > 0
    HFO_rate = numel(slow_cHFOs) / quiet_time_spent;
else
    HFO_rate = NaN;
end

save(swr_file, 'lfp_bp', 'HFO_rate', 'all_cHFOs', 'slow_cHFOs', ...
    'fast_cHFOs', 'quiet_time_spent', 'fast_time_spent', ...
    'quiet_percent', 'fast_percent', 'time_spent', 'detection_threshold');

if save_hfo_rate_file
    save_hfo_rate(processed_folder, sleep_label, HFO_rate);
end

summary.status = 'processed';
summary.slow_event_count = numel(slow_cHFOs);
summary.all_event_count = numel(all_cHFOs);
summary.fast_event_count = numel(fast_cHFOs);
summary.HFO_rate = HFO_rate;
summary.quiet_time_spent = quiet_time_spent;
summary.fast_time_spent = fast_time_spent;
summary.time_spent = time_spent;
end

function indata_idx = resolve_sleep_indata_index(sleep_label, fallback_idx)

tokens = regexp(lower(char(string(sleep_label))), '^s(\d+)$', 'tokens', 'once');
if isempty(tokens)
    indata_idx = fallback_idx;
else
    indata_idx = str2double(tokens{1});
end
end

function save_hfo_rate(processed_folder, sleep_label, HFO_rate)

hfo_rate_file = fullfile(processed_folder, 'HFOrate.mat');
switch lower(sleep_label)
    case 's1'
        i_HFOrate_s1 = [];
        i_HFOrate_s1(2) = HFO_rate;
        save(hfo_rate_file, 'i_HFOrate_s1');
    case 's2'
        i_HFOrate_s2 = [];
        i_HFOrate_s2(2) = HFO_rate;
        save(hfo_rate_file, 'i_HFOrate_s2');
    otherwise
        save(hfo_rate_file, 'HFO_rate');
end
end

function lfp = load_sleep_lfp(main_dir, sleep_label, channel)

eeg_file = fullfile(main_dir, sleep_label, strcat('CSC', num2str(channel), '.ncs'));
fprintf('Reading: %s\n', eeg_file);
[eeg, sample_frequency] = readCRTsd(eeg_file);
lfp.samp = Data(eeg) * -1;
lfp.samp = downsample(lfp.samp, 16);
lfp.ts = Range(eeg);
lfp.ts = downsample(lfp.ts, 16);
lfp.sampFreq = sample_frequency / 16;
end

function [eeg_ts, eeg_raw] = align_lfp_with_indata(tracking_time, eeg_ts, eeg_raw)

if eeg_ts(end) / tracking_time(end) > 10
    eeg_ts = eeg_ts * 1e-4;
end

[~, lfp_start_idx] = min(abs(tracking_time(1) - eeg_ts));
[~, lfp_end_idx] = min(abs(tracking_time(end) - eeg_ts));

eeg_raw = eeg_raw(lfp_start_idx:lfp_end_idx);
eeg_ts = eeg_ts(lfp_start_idx:lfp_end_idx);
end

function values = resize_to_lfp_samples(values, sample_count)

values = double(values(:));
if numel(values) == sample_count
    return
end
if numel(values) == 1
    values = repmat(values, sample_count, 1);
    return
end

source_x = linspace(1, sample_count, numel(values));
target_x = (1:sample_count)';
values = interp1(source_x(:), values(:), target_x, 'linear', 'extrap');
end

function detection_threshold = resolve_detection_threshold(custom_settings)

detection_threshold = struct( ...
    'numSD', 3.5, ...
    'numHighPwrCycles', 4, ...
    'peakSD', 3);

if isfield(custom_settings, 'detectionThreshold') && isstruct(custom_settings.detectionThreshold)
    threshold_override = custom_settings.detectionThreshold;
    fields = fieldnames(threshold_override);
    for field_idx = 1:numel(fields)
        detection_threshold.(fields{field_idx}) = threshold_override.(fields{field_idx});
    end
end

detection_threshold.numSD = double(get_option_value(custom_settings, 'numSD', detection_threshold.numSD));
detection_threshold.numHighPwrCycles = double(get_option_value(custom_settings, 'numHighPwrCycles', detection_threshold.numHighPwrCycles));
detection_threshold.peakSD = double(get_option_value(custom_settings, 'peakSD', detection_threshold.peakSD));
end

function swr_file_name = build_swr_event_filename(detection_threshold)

swr_file_name = sprintf('_allE_numSD%s_HighPwrCycles%s.mat', ...
    num2str(detection_threshold.numSD), num2str(detection_threshold.numHighPwrCycles));
end

function status = summarize_session_status(sleep_summaries)

if any(strcmp({sleep_summaries.status}, 'processed'))
    status = 'processed';
elseif all(strcmp({sleep_summaries.status}, 'skipped_existing_swr_file'))
    status = 'skipped_existing_swr_files';
else
    status = 'skipped';
end
end

function summary = empty_detection_summary()

summary = struct( ...
    'session_index', NaN, ...
    'status', '', ...
    'sleep', [] ...
    );
end

function summary = empty_sleep_detection_summary()

summary = struct( ...
    'sleep_label', '', ...
    'status', '', ...
    'error_message', '', ...
    'swr_file', '', ...
    'slow_event_count', 0, ...
    'all_event_count', 0, ...
    'fast_event_count', 0, ...
    'HFO_rate', NaN, ...
    'quiet_time_spent', NaN, ...
    'fast_time_spent', NaN, ...
    'time_spent', NaN ...
    );
end

function values = normalize_text_cell(values)

if ischar(values) || isstring(values)
    values = cellstr(string(values));
end
end

function config = load_oscillation_config(config_path)

if nargin < 1 || isempty(config_path)
    config_path = fullfile(fileparts(mfilename('fullpath')), '..', 'classification', 'classification_config.json');
end

config = struct('mclustPath', '', 'oscillationAnalysisPath', '');
if exist(config_path, 'file') ~= 2
    return
end

try
    decoded = jsondecode(fileread(config_path));
catch ME
    warning('Could not read config %s: %s', config_path, ME.message);
    return
end

if isfield(decoded, 'mclustPath')
    config.mclustPath = char(string(decoded.mclustPath));
end
if isfield(decoded, 'oscillationAnalysisPath')
    config.oscillationAnalysisPath = char(string(decoded.oscillationAnalysisPath));
end
end

function add_dependency_paths(config, custom_settings)

matlab_code_folder = fullfile(fileparts(mfilename('fullpath')), '..');
if exist(matlab_code_folder, 'dir') == 7
    addpath(genpath(matlab_code_folder));
end

mclust_path = char(string(get_option_value(custom_settings, 'mclustPath', config.mclustPath)));
oscillation_analysis_path = char(string(get_option_value(custom_settings, 'oscillationAnalysisPath', config.oscillationAnalysisPath)));
add_dependency_path(mclust_path, 'MClust');
add_dependency_path(oscillation_analysis_path, 'oscillation analysis dependency');

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

function value = get_option_value(settings_struct, field_name, default_value)

value = default_value;
if nargin < 1 || isempty(settings_struct) || ~isstruct(settings_struct)
    return
end

if isfield(settings_struct, field_name) && ~isempty(settings_struct.(field_name))
    value = settings_struct.(field_name);
end
end
