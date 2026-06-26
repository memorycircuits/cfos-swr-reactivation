function result = compute_swr_phase_analysis(custom_settings)
% Compute per-cell SWR spike phase metrics.

if nargin < 1 || isempty(custom_settings)
    custom_settings = struct();
end

config = load_oscillation_config(get_option_value(custom_settings, 'configPath', ''));
add_dependency_paths(config, custom_settings);

all_cells_path = char(string(get_option_value(custom_settings, 'allCellsPath', '')));
if isempty(all_cells_path)
    all_cells_path = resolve_all_cells_path();
end

session_info_path = char(string(get_option_value(custom_settings, 'sessionInfoPath', '')));
if isempty(session_info_path)
    session_info_path = resolve_session_info_path();
end

loaded_all_cells = load(all_cells_path, 'All_Cells_combined');
if ~isfield(loaded_all_cells, 'All_Cells_combined')
    error('All_Cells_combined was not found in %s.', all_cells_path);
end
All_Cells_combined = loaded_all_cells.All_Cells_combined;

loaded_session_info = load(session_info_path, 'sessInfo');
if ~isfield(loaded_session_info, 'sessInfo')
    error('sessInfo was not found in %s.', session_info_path);
end
sessInfo = loaded_session_info.sessInfo;

session_indices = get_option_value(custom_settings, 'sessionIndices', []);
if isempty(session_indices)
    session_indices = 1:min(numel(sessInfo), numel(All_Cells_combined));
end

sleep_folders = normalize_text_cell(get_option_value(custom_settings, 'sleepFolders', {'s1', 's2'}));
swr_file_name = char(string(get_option_value(custom_settings, 'swrFileName', '_allE_numSD3.5_HighPwrCycles4.mat')));
phase_file_name = char(string(get_option_value(custom_settings, 'phaseFileName', 'Instantaneous_phase.mat')));
save_phase_file = logical(get_option_value(custom_settings, 'savePhaseFile', true));
save_updated_all_cells = logical(get_option_value(custom_settings, 'saveUpdatedAllCells', true));
phase_downsampling_factor = double(get_option_value(custom_settings, 'phaseDownsamplingFactor', 4));

session_summaries = repmat(empty_session_summary(), numel(session_indices), 1);

for session_counter = 1:numel(session_indices)
    session_index = session_indices(session_counter);
    session_summaries(session_counter).session_index = session_index;

    if session_index < 1 || session_index > numel(sessInfo) || session_index > numel(All_Cells_combined)
        session_summaries(session_counter).status = 'skipped_invalid_session_index';
        warning('Skipping invalid session index: %d', session_index);
        continue
    end

    main_dir = normalize_text_value(sessInfo(session_index).mainDir);
    tt_list_path = fullfile(main_dir, normalize_text_value(sessInfo(session_index).tList));
    if exist(tt_list_path, 'file') ~= 2
        session_summaries(session_counter).status = 'skipped_missing_tetrode_list';
        warning('Session %d: missing tList file %s.', session_index, tt_list_path);
        continue
    end

    tt_files = read_tt_list(tt_list_path);
    num_cells = max(numel(tt_files), max_num_cells_in_session(All_Cells_combined(session_index)));
    if num_cells == 0
        session_summaries(session_counter).status = 'skipped_no_cells';
        continue
    end

    sleep_summaries = repmat(empty_sleep_summary(), numel(sleep_folders), 1);
    for sleep_idx = 1:numel(sleep_folders)
        sleep_label = char(string(sleep_folders{sleep_idx}));
        sleep_summaries(sleep_idx) = compute_session_sleep_swr_phase( ...
            sessInfo(session_index), session_index, sleep_label, sleep_idx, ...
            tt_files, num_cells, swr_file_name, phase_file_name, ...
            phase_downsampling_factor, save_phase_file);

        if strcmp(sleep_summaries(sleep_idx).status, 'processed')
            prefix = upper(sleep_label);
            All_Cells_combined = write_phase_metrics_to_all_cells( ...
                All_Cells_combined, session_index, prefix, sleep_summaries(sleep_idx).metrics);
        end
    end

    session_summaries(session_counter).sleep = sleep_summaries;
    session_summaries(session_counter).status = summarize_session_status(sleep_summaries);
end

if save_updated_all_cells
    save(all_cells_path, 'All_Cells_combined', '-append');
end

result = struct( ...
    'allCellsPath', all_cells_path, ...
    'sessionInfoPath', session_info_path, ...
    'swrFileName', swr_file_name, ...
    'phaseFileName', phase_file_name, ...
    'sessions', session_summaries);

fprintf('Computed SWR phase metrics for %d sessions.\n', nnz(strcmp({session_summaries.status}, 'processed')));
end


function summary = compute_session_sleep_swr_phase(sess_entry, session_index, sleep_label, sleep_idx, tt_files, num_cells, swr_file_name, phase_file_name, phase_downsampling_factor, save_phase_file)

summary = empty_sleep_summary();
summary.sleep_label = sleep_label;
summary.swr_file = fullfile(sess_entry.mainDir, sleep_label, 'processedData', swr_file_name);
summary.phase_file = fullfile(sess_entry.mainDir, sleep_label, 'processedData', phase_file_name);

if exist(summary.swr_file, 'file') ~= 2
    summary.status = 'skipped_missing_swr_file';
    warning('Missing SWR event file: %s', summary.swr_file);
    return
end

loaded_swr = load(summary.swr_file, 'slow_cHFOs', 'all_cHFOs');
if ~isfield(loaded_swr, 'slow_cHFOs')
    summary.status = 'skipped_missing_slow_cHFOs';
    warning('SWR file is missing slow_cHFOs: %s', summary.swr_file);
    return
end
slow_cHFOs = loaded_swr.slow_cHFOs;

try
    spike_data = readSpikeDataOnly(fullfile(sess_entry.mainDir, sleep_label), tt_files);
    spike_times = fixSpikes(spike_data);

    indata_file = fullfile(sess_entry.mainDir, 'processedData', 'indataS.mat');
    loaded_indata = load(indata_file, 'indata');
    indata_idx = resolve_sleep_indata_index(sleep_label, sleep_idx);
    if ~isfield(loaded_indata, 'indata') || isnan(indata_idx) || indata_idx < 1 || indata_idx > numel(loaded_indata.indata)
        error('Could not load %s indata index %d from %s.', sleep_label, indata_idx, indata_file);
    end
    indata = loaded_indata.indata(indata_idx);

    lfp = load_sleep_lfp(sess_entry.mainDir, sleep_label, sess_entry.cellLayerChann, phase_downsampling_factor);
    [lfp.ts, lfp.samp] = align_lfp_with_indata(indata.t, lfp.ts, lfp.samp);
    [lfp_zscore, inst_phase] = compute_ripple_phase_trace(lfp);

    metrics = compute_sleep_phase_metrics(lfp, lfp_zscore, inst_phase, slow_cHFOs, spike_times, num_cells);
catch ME
    summary.status = 'failed';
    summary.error_message = ME.message;
    warning('SWR phase analysis failed for %s %s: %s', sess_entry.mainDir, sleep_label, ME.message);
    return
end

phase_metadata = build_phase_metadata(sess_entry, session_index, sleep_label, tt_files, num_cells, numel(slow_cHFOs), phase_downsampling_factor, swr_file_name);
metrics.phase_metadata = phase_metadata;

if save_phase_file
    save_swr_phase_file(summary.phase_file, metrics);
end

summary.status = 'processed';
summary.num_cells = num_cells;
summary.num_swrs = numel(slow_cHFOs);
summary.metrics = metrics;
end


function [lfp_zscore, inst_phase] = compute_ripple_phase_trace(lfp)

Hd = detectHFOs.blanco_bp(lfp.sampFreq);
filter_numerator = Hd.Numerator;
lfp_bp = filtfilt(filter_numerator, 1, lfp.samp);
lfp_zscore = zscore_vector(lfp_bp);
analytic_signal = hilbert(lfp_zscore);
inst_phase = angle(analytic_signal);
end


function metrics = compute_sleep_phase_metrics(lfp, lfp_zscore, inst_phase, slow_cHFOs, spike_times, num_cells)

num_swrs = numel(slow_cHFOs);
cell_inst_phase_in_swr = cell(num_swrs, num_cells);
cell_inst_phase_aligned_to_trough = cell(num_swrs, num_cells);
cell_center_aligned_spikes_phase = cell(num_swrs, num_cells);

for swr_idx = 1:num_swrs
    [swr_start_idx, swr_stop_idx] = swr_time_to_lfp_indices(lfp.ts, slow_cHFOs(swr_idx));
    if swr_stop_idx <= swr_start_idx
        continue
    end

    swr_signal = lfp_zscore(swr_start_idx:swr_stop_idx);
    swr_timestamps = lfp.ts(swr_start_idx:swr_stop_idx);
    [mean_trough_phase, swr_phase_trace, has_center_phase] = compute_swr_alignment_phase(swr_signal, inst_phase, swr_start_idx);

    swr_start = swr_timestamps(1);
    swr_stop = swr_timestamps(end);
    for cell_idx = 1:min(numel(spike_times), num_cells)
        cell_spikes = spike_times{cell_idx};
        if isempty(cell_spikes)
            continue
        end

        cell_spikes = double(cell_spikes(:));
        spikes_in_swr = cell_spikes(cell_spikes >= swr_start & cell_spikes <= swr_stop);
        if isempty(spikes_in_swr)
            continue
        end

        lfp_idx_rel = nearest_indices(swr_timestamps, spikes_in_swr);
        lfp_idx_abs = swr_start_idx + lfp_idx_rel - 1;
        spikes_phase = inst_phase(lfp_idx_abs);
        cell_inst_phase_in_swr{swr_idx, cell_idx} = spikes_phase(:)';

        if isfinite(mean_trough_phase)
            cell_inst_phase_aligned_to_trough{swr_idx, cell_idx} = wrap_to_pi(spikes_phase - mean_trough_phase);
            cell_inst_phase_aligned_to_trough{swr_idx, cell_idx} = cell_inst_phase_aligned_to_trough{swr_idx, cell_idx}(:)';
        end
        if has_center_phase
            cell_center_aligned_spikes_phase{swr_idx, cell_idx} = swr_phase_trace(lfp_idx_rel);
            cell_center_aligned_spikes_phase{swr_idx, cell_idx} = cell_center_aligned_spikes_phase{swr_idx, cell_idx}(:)';
        end
    end
end

[cells_phase_mean_angle, cells_phase_R, cells_phase_p_val, cells_phase_z] = ...
    compute_cell_phase_statistics(cell_inst_phase_in_swr, num_cells);
[cells_phase_mean_angle_troughAligned, cells_phase_R_troughAligned, ...
    cells_phase_p_val_troughAligned, cells_phase_z_troughAligned] = ...
    compute_cell_phase_statistics(cell_inst_phase_aligned_to_trough, num_cells);

metrics = struct( ...
    'cells_phase_mean_angle', cells_phase_mean_angle, ...
    'cells_phase_R', cells_phase_R, ...
    'cells_phase_p_val', cells_phase_p_val, ...
    'cells_phase_z', cells_phase_z, ...
    'cells_phase_mean_angle_troughAligned', cells_phase_mean_angle_troughAligned, ...
    'cells_phase_R_troughAligned', cells_phase_R_troughAligned, ...
    'cells_phase_p_val_troughAligned', cells_phase_p_val_troughAligned, ...
    'cells_phase_z_troughAligned', cells_phase_z_troughAligned, ...
    'cell_inst_phase_in_swr', {cell_inst_phase_in_swr}, ...
    'cell_inst_phase_aligned_to_trough', {cell_inst_phase_aligned_to_trough}, ...
    'cell_center_aligned_spikes_phase', {cell_center_aligned_spikes_phase});
end


function [mean_trough_phase, swr_phase_trace, has_center_phase] = compute_swr_alignment_phase(swr_signal, inst_phase, swr_start_idx)

[~, trough_locs_rel] = findpeaks(-swr_signal);
if isempty(trough_locs_rel)
    mean_trough_phase = NaN;
else
    trough_locs_abs = swr_start_idx + trough_locs_rel - 1;
    mean_trough_phase = circular_mean(inst_phase(trough_locs_abs));
end

[~, peak_locs_rel] = findpeaks(swr_signal);
if isempty(peak_locs_rel)
    swr_phase_trace = NaN(size(swr_signal));
    has_center_phase = false;
    return
end

swr_center_idx = (numel(swr_signal) + 1) / 2;
[~, center_peak_idx] = min(abs(peak_locs_rel - swr_center_idx));
center_peak_loc_rel = peak_locs_rel(center_peak_idx);

raw_inst_phase = angle(hilbert(swr_signal));
unwrapped_phase = unwrap(raw_inst_phase);
swr_phase_trace = unwrapped_phase - unwrapped_phase(center_peak_loc_rel);
has_center_phase = true;
end


function [mean_angle, phase_R, p_val, z_val] = compute_cell_phase_statistics(phase_cell_array, num_cells)

mean_angle = NaN(num_cells, 1);
phase_R = NaN(num_cells, 1);
p_val = NaN(num_cells, 1);
z_val = NaN(num_cells, 1);

for cell_idx = 1:num_cells
    phases = concatenate_phase_cells(phase_cell_array(:, cell_idx));
    if isempty(phases)
        continue
    end

    phases = phases(:);
    phases = phases(isfinite(phases));
    if isempty(phases)
        continue
    end

    mean_angle(cell_idx) = circular_mean(phases);
    phase_R(cell_idx) = circular_resultant_length(phases);
    [p_val(cell_idx), z_val(cell_idx)] = circular_rayleigh_test(phases);
end
end


function phases = concatenate_phase_cells(phase_cells)

phase_cells = phase_cells(~cellfun('isempty', phase_cells));
if isempty(phase_cells)
    phases = [];
else
    phases = horzcat(phase_cells{:});
end
end


function all_cells = write_phase_metrics_to_all_cells(all_cells, session_index, prefix, metrics)

all_cells(session_index).(sprintf('%s_cells_phase_mean_angle', prefix)) = metrics.cells_phase_mean_angle;
all_cells(session_index).(sprintf('%s_cells_phase_R', prefix)) = metrics.cells_phase_R;
all_cells(session_index).(sprintf('%s_cells_phase_p_val', prefix)) = metrics.cells_phase_p_val;
all_cells(session_index).(sprintf('%s_cells_phase_z', prefix)) = metrics.cells_phase_z;

all_cells(session_index).(sprintf('%s_cells_phase_mean_angle_troughAligned', prefix)) = metrics.cells_phase_mean_angle_troughAligned;
all_cells(session_index).(sprintf('%s_cells_phase_R_troughAligned', prefix)) = metrics.cells_phase_R_troughAligned;
all_cells(session_index).(sprintf('%s_cells_phase_p_val_troughAligned', prefix)) = metrics.cells_phase_p_val_troughAligned;
all_cells(session_index).(sprintf('%s_cells_phase_z_troughAligned', prefix)) = metrics.cells_phase_z_troughAligned;
end


function save_swr_phase_file(phase_file, metrics)

phase_folder = fileparts(phase_file);
if exist(phase_folder, 'dir') ~= 7
    mkdir(phase_folder);
end

cells_phase_mean_angle = metrics.cells_phase_mean_angle;
cells_phase_R = metrics.cells_phase_R;
cells_phase_p_val = metrics.cells_phase_p_val;
cells_phase_z = metrics.cells_phase_z;
cells_phase_mean_angle_troughAligned = metrics.cells_phase_mean_angle_troughAligned;
cells_phase_R_troughAligned = metrics.cells_phase_R_troughAligned;
cells_phase_p_val_troughAligned = metrics.cells_phase_p_val_troughAligned;
cells_phase_z_troughAligned = metrics.cells_phase_z_troughAligned;
cell_inst_phase_in_swr = metrics.cell_inst_phase_in_swr;
cell_inst_phase_aligned_to_trough = metrics.cell_inst_phase_aligned_to_trough;
cell_center_aligned_spikes_phase = metrics.cell_center_aligned_spikes_phase;
phase_metadata = metrics.phase_metadata;

save(phase_file, 'cells_phase_mean_angle', 'cells_phase_R', 'cells_phase_p_val', 'cells_phase_z', ...
    'cells_phase_mean_angle_troughAligned', 'cells_phase_R_troughAligned', ...
    'cells_phase_p_val_troughAligned', 'cells_phase_z_troughAligned', ...
    'cell_inst_phase_in_swr', 'cell_inst_phase_aligned_to_trough', ...
    'cell_center_aligned_spikes_phase', 'phase_metadata');
end


function lfp = load_sleep_lfp(main_dir, sleep_label, channel, downsampling_factor)

if iscell(channel)
    channel = channel{1};
end
if ischar(channel) || isstring(channel)
    channel = str2double(char(string(channel)));
end
channel = double(channel(1));

eeg_file = fullfile(main_dir, sleep_label, strcat('CSC', num2str(channel), '.ncs'));
fprintf('Reading: %s\n', eeg_file);
[eeg, sample_frequency] = readCRTsd(eeg_file);
lfp.samp = Data(eeg) * -1;
lfp.samp = downsample(lfp.samp(:), downsampling_factor);
lfp.ts = Range(eeg);
lfp.ts = downsample(lfp.ts(:), downsampling_factor);
lfp.sampFreq = sample_frequency / downsampling_factor;
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


function [start_idx, stop_idx] = swr_time_to_lfp_indices(lfp_timestamps, swr_event)

[~, start_idx] = min(abs(lfp_timestamps - swr_event.start_ts));
[~, stop_idx] = min(abs(lfp_timestamps - swr_event.stop_ts));
start_idx = max(1, min(start_idx, numel(lfp_timestamps)));
stop_idx = max(1, min(stop_idx, numel(lfp_timestamps)));
end


function idx = nearest_indices(reference_values, query_values)

idx = zeros(numel(query_values), 1);
for query_idx = 1:numel(query_values)
    [~, idx(query_idx)] = min(abs(reference_values - query_values(query_idx)));
end
end


function lfp_zscore = zscore_vector(values)

values = double(values(:));
mu = mean(values(isfinite(values)));
sigma = std(values(isfinite(values)));
if ~isfinite(sigma) || sigma == 0
    sigma = 1;
end
lfp_zscore = (values - mu) ./ sigma;
end


function mean_angle = circular_mean(phases)

phases = phases(:);
mean_angle = atan2(mean(sin(phases)), mean(cos(phases)));
end


function resultant_length = circular_resultant_length(phases)

phases = phases(:);
resultant_length = sqrt(sum(sin(phases)).^2 + sum(cos(phases)).^2) / numel(phases);
end


function [p_val, z_val] = circular_rayleigh_test(phases)

phases = phases(:);
n = numel(phases);
r = circular_resultant_length(phases);
z_val = n * r.^2;
R = n * r;
p_val = exp(sqrt(1 + 4*n + 4*(n^2 - R^2)) - (1 + 2*n));

if n < 50
    if z_val < 10
        p_val = exp(-z_val) * (1 + (2*z_val - z_val^2) / (4*n) - ...
            (24*z_val - 132*z_val^2 + 76*z_val^3 - 9*z_val^4) / (288*n^2));
    else
        p_val = exp(-z_val) * (1 - (1 + 2*z_val) / (4*n) + ...
            (1 + 2*z_val + 4*z_val^2) / (8*n^2));
    end
end

p_val = min(max(p_val, 0), 1);
end


function wrapped = wrap_to_pi(phases)

wrapped = mod(phases + pi, 2*pi) - pi;
end


function tt_files = read_tt_list(tt_list_path)

fid = fopen(tt_list_path);
if fid < 0
    error('Could not open tList file: %s', tt_list_path);
end

tt_files = {};
cleanup = onCleanup(@() fclose(fid));
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
delete(cleanup);
end


function indata_idx = resolve_sleep_indata_index(sleep_label, fallback_idx)

tokens = regexp(lower(char(string(sleep_label))), '^s(\d+)$', 'tokens', 'once');
if isempty(tokens)
    indata_idx = fallback_idx;
else
    indata_idx = str2double(tokens{1});
end
end


function metadata = build_phase_metadata(sess_entry, session_index, sleep_label, tt_files, num_cells, num_swrs, phase_downsampling_factor, swr_file_name)

metadata = struct();
metadata.session_index = session_index;
metadata.animal = get_struct_field(sess_entry, 'animal');
metadata.day = get_struct_field(sess_entry, 'day');
metadata.sleep_folder = sleep_label;
metadata.num_cells = num_cells;
metadata.num_swrs = num_swrs;
metadata.t_list = tt_files(:);
metadata.t_list_file = get_struct_field(sess_entry, 'tList');
metadata.cell_layer_channel = get_struct_field(sess_entry, 'cellLayerChann');
metadata.swr_file_name = swr_file_name;
metadata.phase_downsampling_factor = phase_downsampling_factor;
metadata.phase_convention = 'angle(hilbert(zscored ripple-band LFP)); trough-aligned phase subtracts circular mean trough phase per SWR';
metadata.source_script = mfilename;
metadata.generated_on = char(datetime('now', 'Format', 'yyyy-MM-dd HH:mm:ss'));
end


function value = get_struct_field(input_struct, field_name)

if isfield(input_struct, field_name)
    value = input_struct.(field_name);
else
    value = [];
end
end


function num_cells = max_num_cells_in_session(session_cells)

num_cells = 0;
candidate_fields = {'S1_cells_phase_mean_angle', 'S2_cells_phase_mean_angle', ...
    'S1_PSP', 'S2_PSP', 'optotagged', 'final_classification_numeric', ...
    'GMM_based_classification_days', 'thetaPhase_pref', 'of_avg_fir_rate'};

for field_idx = 1:numel(candidate_fields)
    field_name = candidate_fields{field_idx};
    if isfield(session_cells, field_name) && ~isempty(session_cells.(field_name))
        num_cells = max(num_cells, numel(session_cells.(field_name)));
    end
end
end


function status = summarize_session_status(sleep_summaries)

if any(strcmp({sleep_summaries.status}, 'processed'))
    status = 'processed';
elseif all(strcmp({sleep_summaries.status}, 'skipped_missing_swr_file'))
    status = 'skipped_missing_swr_files';
else
    status = 'skipped';
end
end


function summary = empty_session_summary()

summary = struct( ...
    'session_index', NaN, ...
    'status', '', ...
    'sleep', [] ...
    );
end


function summary = empty_sleep_summary()

summary = struct( ...
    'sleep_label', '', ...
    'status', '', ...
    'error_message', '', ...
    'swr_file', '', ...
    'phase_file', '', ...
    'num_cells', 0, ...
    'num_swrs', 0, ...
    'metrics', [] ...
    );
end


function text_value = normalize_text_value(value)

if iscell(value)
    value = value{1};
end
if isstring(value)
    text_value = char(value(1));
elseif ischar(value)
    text_value = value;
else
    text_value = char(string(value));
end
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

config = struct('mclustPath', '', 'oscillationAnalysisPath', '', 'swrPhaseAnalysisPath', '');
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
if isfield(decoded, 'swrPhaseAnalysisPath')
    config.swrPhaseAnalysisPath = char(string(decoded.swrPhaseAnalysisPath));
end
end


function add_dependency_paths(config, custom_settings)

matlab_code_folder = fullfile(fileparts(mfilename('fullpath')), '..');
if exist(matlab_code_folder, 'dir') == 7
    addpath(genpath(matlab_code_folder));
end

mclust_path = char(string(get_option_value(custom_settings, 'mclustPath', config.mclustPath)));
oscillation_analysis_path = char(string(get_option_value(custom_settings, 'oscillationAnalysisPath', config.oscillationAnalysisPath)));
swr_phase_analysis_path = char(string(get_option_value(custom_settings, 'swrPhaseAnalysisPath', config.swrPhaseAnalysisPath)));

add_dependency_path(mclust_path, 'MClust');
add_dependency_path(oscillation_analysis_path, 'oscillation analysis dependency');
add_dependency_path(swr_phase_analysis_path, 'SWR phase analysis dependency');

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
