function result = compute_swr_cell_metrics(custom_settings)
% Compute per-cell SWR metrics stored in All_Cells_combined.

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
swr_file_name = char(string(get_option_value(custom_settings, 'swrFileName', '_allE_numSD3.5_HighPwrCycles4.mat')));
save_updated_all_cells = logical(get_option_value(custom_settings, 'saveUpdatedAllCells', true));
save_swr_data = logical(get_option_value(custom_settings, 'saveSWRData', true));
update_script_name = 'compute_swr_cell_metrics.m';

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

session_summaries = repmat(empty_swr_summary(), numel(session_indices), 1);

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

    processed_sleeps = {};
    for sleep_idx = 1:numel(sleep_folders)
        sleep_label = sleep_folders{sleep_idx};
        metrics = compute_session_sleep_swr_metrics( ...
            sessInfo(session_index), main_dir, sleep_label, sleep_idx, tt_files, swr_file_name);

        if ~strcmp(metrics.status, 'processed')
            session_summaries(idx).sleep_statuses{sleep_idx} = metrics.status;
            continue
        end

        prefix = upper(sleep_label);
        [All_Cells_combined, written_fields] = write_swr_metrics_to_all_cells( ...
            All_Cells_combined, session_index, prefix, metrics);
        last_updated = mark_fields_updated(last_updated, session_index, written_fields, update_script_name);

        if save_swr_data
            save_swr_data_file(metrics.swrDataFile, metrics);
        end

        processed_sleeps{end + 1} = sleep_label; %#ok<AGROW>
        session_summaries(idx).sleep_statuses{sleep_idx} = 'processed';
        session_summaries(idx).swr_counts(sleep_idx) = metrics.SWR_number;
    end

    if isempty(processed_sleeps)
        session_summaries(idx).status = 'skipped_no_processed_sleep_sessions';
    else
        session_summaries(idx).status = 'processed';
        session_summaries(idx).processed_sleeps = processed_sleeps;
        session_summaries(idx).cell_count = numel(tt_files);
    end
end

if save_updated_all_cells
    save_updated_all_cells_file(all_cells_path, All_Cells_combined, last_updated);
end

result = struct( ...
    'allCellsPath', all_cells_path, ...
    'sessionInfoPath', session_info_path, ...
    'swrFileName', swr_file_name, ...
    'sessions', session_summaries);

fprintf('Computed SWR cell metrics for %d sessions.\n', ...
    nnz(strcmp({session_summaries.status}, 'processed')));
end

function metrics = compute_session_sleep_swr_metrics(sess_entry, main_dir, sleep_label, sleep_idx, tt_files, swr_file_name)

metrics = empty_metrics();
metrics.sleepLabel = sleep_label;
metrics.swrDataFile = fullfile(main_dir, sleep_label, 'processedData', 'SWR_data.mat');

[spike_times, loaded_ok] = read_subsession_spikes(fullfile(main_dir, sleep_label), tt_files);
if ~loaded_ok
    metrics.status = 'skipped_unreadable_spikes';
    return
end

indata_file = fullfile(main_dir, 'processedData', 'indataS.mat');
if exist(indata_file, 'file') ~= 2
    metrics.status = 'skipped_missing_indataS';
    warning('Missing indataS file: %s', indata_file);
    return
end
loaded_indata = load(indata_file, 'indata');
if ~isfield(loaded_indata, 'indata') || sleep_idx > numel(loaded_indata.indata)
    metrics.status = 'skipped_invalid_indataS';
    warning('Could not load sleep %d indata from %s.', sleep_idx, indata_file);
    return
end
indata = loaded_indata.indata(sleep_idx);

swr_file = fullfile(main_dir, sleep_label, 'processedData', swr_file_name);
if exist(swr_file, 'file') ~= 2
    metrics.status = 'skipped_missing_swr_file';
    warning('Missing SWR detection file: %s', swr_file);
    return
end
loaded_swr = load(swr_file, 'slow_cHFOs');
if ~isfield(loaded_swr, 'slow_cHFOs')
    metrics.status = 'skipped_missing_slow_cHFOs';
    warning('SWR file did not contain slow_cHFOs: %s', swr_file);
    return
end
slow_cHFOs = loaded_swr.slow_cHFOs;

lfp = load_sleep_lfp(main_dir, sleep_label, sess_entry.cellLayerChann);
[lfp.ts, lfp.samp] = align_lfp_with_indata(indata.t, lfp.ts, lfp.samp);
session_duration = lfp.ts(end) - lfp.ts(1);

num_cells = numel(tt_files);
num_swrs = numel(slow_cHFOs);

SWRspikes = zeros(num_swrs, num_cells);
pre_SWRspikes = zeros(num_swrs, num_cells);
SWRrealspikes = zeros(num_swrs, num_cells);
pre_SWR_real_spikes = zeros(num_swrs, num_cells);
SpPR_SWR = NaN(num_swrs, num_cells);
SSMI_SWR = NaN(num_swrs, num_cells);
time_of_spikes = cell(num_swrs, num_cells);
mean_time_of_spikes = zeros(num_swrs, num_cells);
median_time_of_spikes = zeros(num_swrs, num_cells);
first_time_of_spikes = zeros(num_swrs, num_cells);
SWR_spiking_activity = cell(num_cells, 1);
SWR_activation_durations = cell(num_cells, 1);
SWR_modulation = NaN(1, num_cells);
unfiltered_SWR_modulation = NaN(1, num_cells);
SWR_modulation_real_duration = NaN(1, num_cells);
swr_duration = zeros(num_swrs, 1);
peaks_amp_sd = zeros(num_swrs, 1);

all_starts = [slow_cHFOs.start_ts];
all_stops = [slow_cHFOs.stop_ts];

for swr_idx = 1:num_swrs
    swr_start = slow_cHFOs(swr_idx).start_ts;
    swr_stop = slow_cHFOs(swr_idx).stop_ts;
    swr_duration(swr_idx) = swr_stop - swr_start;
    if isfield(slow_cHFOs, 'peaks_sd')
        peaks_amp_sd(swr_idx) = slow_cHFOs(swr_idx).peaks_sd;
    else
        peaks_amp_sd(swr_idx) = NaN;
    end

    zero_point = swr_start + (swr_stop - swr_start) / 2;
    mod_start = zero_point - 0.05;
    mod_stop = zero_point + 0.05;
    swr_duration_timing_start = swr_start - 2 * swr_duration(swr_idx);
    swr_duration_timing_stop = swr_stop + 2 * swr_duration(swr_idx);
    pre_swr_real_start = swr_start - 2 * swr_duration(swr_idx);
    pre_swr_real_stop = swr_start - swr_duration(swr_idx);
    pre_mod_start = zero_point - 0.150;
    pre_mod_stop = zero_point - 0.05;

    shift_step = 0.025;
    max_shift = 1;
    total_shift = 0;
    while true
        overlap = window_overlaps_events(pre_mod_start, pre_mod_stop, all_starts, all_stops);
        if ~overlap
            break
        end
        pre_mod_start = pre_mod_start - shift_step;
        pre_mod_stop = pre_mod_stop - shift_step;
        total_shift = total_shift + shift_step;
        if total_shift >= max_shift
            warning('Could not resolve SWR overlap within max allowed shift.');
            break
        end
    end

    while true
        overlap = window_overlaps_events(pre_swr_real_start, pre_swr_real_stop, all_starts, all_stops);
        if ~overlap
            break
        end
        pre_swr_real_start = pre_swr_real_start - shift_step;
        pre_swr_real_stop = pre_swr_real_stop - shift_step;
        total_shift = total_shift + shift_step;
        if total_shift >= max_shift
            warning('Could not resolve SWR overlap within max allowed shift.');
            break
        end
    end

    for cell_idx = 1:num_cells
        cell_spikes = spike_times{cell_idx};
        if isempty(cell_spikes)
            continue
        end

        spikes_in_swr = find(mod_start <= cell_spikes & cell_spikes <= mod_stop);
        spikes_in_real_swr = find(swr_start <= cell_spikes & cell_spikes <= swr_stop);
        spikes_in_swr_duration = find(swr_duration_timing_start <= cell_spikes & cell_spikes <= swr_duration_timing_stop);

        if ~isempty(spikes_in_swr_duration)
            time_of_spikes{swr_idx, cell_idx} = ...
                (cell_spikes(spikes_in_swr_duration) - swr_start) ./ swr_duration(swr_idx);
            mean_time_of_spikes(swr_idx, cell_idx) = mean(time_of_spikes{swr_idx, cell_idx});
            median_time_of_spikes(swr_idx, cell_idx) = median(time_of_spikes{swr_idx, cell_idx});
            spike_array = time_of_spikes{swr_idx, cell_idx};
            first_time_of_spikes(swr_idx, cell_idx) = spike_array(1);

            current_spikes = (cell_spikes(spikes_in_swr_duration) - swr_start) ./ swr_duration(swr_idx);
            SWR_spiking_activity{cell_idx} = [SWR_spiking_activity{cell_idx}, reshape(current_spikes, 1, [])];
            SWR_activation_durations{cell_idx, 1} = [SWR_activation_durations{cell_idx, 1}; swr_duration(swr_idx)];
        end

        if ~isempty(spikes_in_real_swr)
            SpPR_SWR(swr_idx, cell_idx) = numel(spikes_in_real_swr);
        end

        spikes_before_swr = find(pre_mod_start <= cell_spikes & cell_spikes <= pre_mod_stop);
        spikes_before_real_swr = find(pre_swr_real_start <= cell_spikes & cell_spikes <= pre_swr_real_stop);

        SWRspikes(swr_idx, cell_idx) = numel(spikes_in_swr);
        pre_SWRspikes(swr_idx, cell_idx) = numel(spikes_before_swr);
        SWRrealspikes(swr_idx, cell_idx) = numel(spikes_in_real_swr);
        pre_SWR_real_spikes(swr_idx, cell_idx) = numel(spikes_before_real_swr);

        epsilon = 0.0000001;
        SSMI_SWR(swr_idx, cell_idx) = ...
            (SWRrealspikes(swr_idx, cell_idx) - pre_SWR_real_spikes(swr_idx, cell_idx)) / ...
            (SWRrealspikes(swr_idx, cell_idx) + pre_SWR_real_spikes(swr_idx, cell_idx) + epsilon);
    end
end

SWR_sums = sum(SWRspikes, 1);
pre_SWR_sums = sum(pre_SWRspikes, 1);
SWR_real_sums = sum(SWRrealspikes, 1);
pre_SWR_real_sums = sum(pre_SWR_real_spikes, 1);

total_cell_spikes = cellfun(@(spikes) size(spikes, 1), spike_times);
if ~isempty(SWR_real_sums)
    outsideSWR_real_sums = total_cell_spikes - SWR_real_sums';
else
    outsideSWR_real_sums = total_cell_spikes;
end

if num_swrs > 0
    SWR_freq = SWR_sums / (num_swrs * 0.1);
    pre_SWR_freq = pre_SWR_sums / (num_swrs * 0.1);
else
    SWR_freq = NaN(1, num_cells);
    pre_SWR_freq = NaN(1, num_cells);
end

if ~isempty(slow_cHFOs)
    totalSWRs_duration = sum_swr_duration_field(slow_cHFOs, swr_duration);
    SWR_real_freq = SWR_real_sums / totalSWRs_duration;
    total_outsideSWRs_duration = session_duration - totalSWRs_duration;
    outsideSWR_freq = outsideSWR_real_sums' / total_outsideSWRs_duration;
    pre_SWR_real_freq = pre_SWR_real_sums / totalSWRs_duration;
else
    SWR_real_freq = NaN(1, num_cells);
    pre_SWR_real_freq = NaN(1, num_cells);
    outsideSWR_freq = NaN(1, num_cells);
end

PSP = NaN(1, num_cells);
longPSP = NaN(1, num_cells);
SFI = NaN(1, num_cells);
SRR = NaN(1, num_cells);
if isempty(SSMI_SWR)
    SSMI = NaN(1, num_cells);
    SpPR = NaN(1, num_cells);
else
    SSMI = mean(SSMI_SWR, 1);
    SpPR = column_nanmean(SpPR_SWR);
end

epsilon = 0.0000001;
for cell_idx = 1:num_cells
    if ~isempty(SWRrealspikes)
        PSP(cell_idx) = (sum(SWRrealspikes(:, cell_idx) > 0) / numel(SWRrealspikes(:, cell_idx))) * 100;
        long_mask = swr_duration(:) > 0.07;
        longPSP(cell_idx) = (sum(SWRrealspikes(long_mask, cell_idx) > 0) / sum(long_mask)) * 100;
    end

    SRR(cell_idx) = log(SWR_real_freq(cell_idx) / (outsideSWR_freq(cell_idx) + epsilon));
    SFI(cell_idx) = SWR_real_freq(cell_idx) - pre_SWR_real_freq(cell_idx);

    if ~isempty(SWRspikes) && (sum(SWRspikes(:, cell_idx) == 0) / numel(SWRspikes(:, cell_idx))) * 100 >= 90
        SWR_modulation(cell_idx) = NaN;
        if SWR_freq(cell_idx) > pre_SWR_freq(cell_idx)
            unfiltered_SWR_modulation(cell_idx) = (SWR_freq(cell_idx) - pre_SWR_freq(cell_idx)) / SWR_freq(cell_idx);
        else
            unfiltered_SWR_modulation(cell_idx) = (SWR_freq(cell_idx) - pre_SWR_freq(cell_idx)) / pre_SWR_freq(cell_idx);
        end
    elseif SWR_freq(cell_idx) > pre_SWR_freq(cell_idx)
        SWR_modulation(cell_idx) = (SWR_freq(cell_idx) - pre_SWR_freq(cell_idx)) / SWR_freq(cell_idx);
        unfiltered_SWR_modulation(cell_idx) = (SWR_freq(cell_idx) - pre_SWR_freq(cell_idx)) / SWR_freq(cell_idx);
    else
        SWR_modulation(cell_idx) = (SWR_freq(cell_idx) - pre_SWR_freq(cell_idx)) / pre_SWR_freq(cell_idx);
        unfiltered_SWR_modulation(cell_idx) = (SWR_freq(cell_idx) - pre_SWR_freq(cell_idx)) / pre_SWR_freq(cell_idx);
    end

    if SWR_real_freq(cell_idx) > pre_SWR_real_freq(cell_idx)
        SWR_modulation_real_duration(cell_idx) = ...
            (SWR_real_freq(cell_idx) - pre_SWR_real_freq(cell_idx)) / SWR_real_freq(cell_idx);
    else
        SWR_modulation_real_duration(cell_idx) = ...
            (SWR_real_freq(cell_idx) - pre_SWR_real_freq(cell_idx)) / pre_SWR_real_freq(cell_idx);
    end
end

SWR_all_binary = SWRspikes >= 1;
SWR_all_sums = sum(SWR_all_binary, 2);
SWR_all_percentage_reactivation = SWR_all_sums / numel(spike_times);

non_zero_mean_spike_time = arrayfun(@(col) ...
    mean(mean_time_of_spikes(mean_time_of_spikes(:, col) ~= 0, col)), ...
    1:size(mean_time_of_spikes, 2));
non_zero_median_spike_time = arrayfun(@(col) ...
    median(median_time_of_spikes(median_time_of_spikes(:, col) ~= 0, col)), ...
    1:size(median_time_of_spikes, 2));
non_zero_mean_first_time_of_spikes = arrayfun(@(col) ...
    mean(first_time_of_spikes(first_time_of_spikes(:, col) ~= 0, col)), ...
    1:size(first_time_of_spikes, 2));

metrics.status = 'processed';
metrics.SWR_number = num_swrs;
metrics.swr_duration_mean = mean(swr_duration);
metrics.SWR_spiking_activity = SWR_spiking_activity;
metrics.SWR_activation_durations = SWR_activation_durations;
metrics.non_zero_mean_first_time_of_spikes = non_zero_mean_first_time_of_spikes;
metrics.non_zero_median_spike_time = non_zero_median_spike_time;
metrics.non_zero_mean_spike_time = non_zero_mean_spike_time;
metrics.time_of_spikes = time_of_spikes;
metrics.SWRspikes = SWRspikes;
metrics.SWRrealspikes = SWRrealspikes;
metrics.pre_SWR_real_spikes = pre_SWR_real_spikes;
metrics.pre_SWRspikes = pre_SWRspikes;
metrics.SWR_sums = SWR_sums;
metrics.SWR_real_sums = SWR_real_sums;
metrics.pre_SWR_real_sums = pre_SWR_real_sums;
metrics.pre_SWR_sums = pre_SWR_sums;
metrics.SWR_modulation = SWR_modulation;
metrics.unfiltered_SWR_modulation = unfiltered_SWR_modulation;
metrics.SWR_modulation_real_duration = SWR_modulation_real_duration;
metrics.swr_duration = swr_duration;
metrics.SWR_all_percentage_reactivation = SWR_all_percentage_reactivation;
metrics.SWR_all_binary = SWR_all_binary;
metrics.peaks_amp_sd = peaks_amp_sd;
metrics.PSP = PSP;
metrics.SFI = SFI;
metrics.SSMI = SSMI;
metrics.longPSP = longPSP;
metrics.SpPR = SpPR;
metrics.SRR = SRR;
end

function [All_Cells_combined, written_fields] = write_swr_metrics_to_all_cells(All_Cells_combined, session_index, prefix, metrics)

field_values = { ...
    sprintf('%s_SWR_modulation', prefix), metrics.SWR_modulation'; ...
    sprintf('%s_unfiltered_SWR_modulation', prefix), metrics.unfiltered_SWR_modulation'; ...
    sprintf('%s_SWR_modulation_real_duration', prefix), metrics.SWR_modulation_real_duration'; ...
    sprintf('%s_SWR_mean_spikeT', prefix), metrics.non_zero_mean_spike_time'; ...
    sprintf('%s_SWR_median_spikeT', prefix), metrics.non_zero_median_spike_time'; ...
    sprintf('%s_PSP', prefix), metrics.PSP'; ...
    sprintf('%s_SFI', prefix), metrics.SFI'; ...
    sprintf('%s_SSMI', prefix), metrics.SSMI'; ...
    sprintf('%s_longPSP', prefix), metrics.longPSP'; ...
    sprintf('%s_SpPR', prefix), metrics.SpPR'; ...
    sprintf('%s_SRR', prefix), metrics.SRR' ...
    };

written_fields = field_values(:, 1)';
for field_idx = 1:size(field_values, 1)
    All_Cells_combined(session_index).(field_values{field_idx, 1}) = field_values{field_idx, 2};
end
end

function save_swr_data_file(swr_data_file, metrics)

swr_duration_mean = metrics.swr_duration_mean;
SWR_spiking_activity = metrics.SWR_spiking_activity;
SWR_activation_durations = metrics.SWR_activation_durations;
SWR_number = metrics.SWR_number;
non_zero_mean_first_time_of_spikes = metrics.non_zero_mean_first_time_of_spikes;
non_zero_median_spike_time = metrics.non_zero_median_spike_time;
non_zero_mean_spike_time = metrics.non_zero_mean_spike_time;
time_of_spikes = metrics.time_of_spikes;
SWRspikes = metrics.SWRspikes;
SWRrealspikes = metrics.SWRrealspikes;
pre_SWR_real_spikes = metrics.pre_SWR_real_spikes;
pre_SWRspikes = metrics.pre_SWRspikes;
SWR_sums = metrics.SWR_sums;
SWR_real_sums = metrics.SWR_real_sums;
pre_SWR_real_sums = metrics.pre_SWR_real_sums;
pre_SWR_sums = metrics.pre_SWR_sums;
SWR_modulation = metrics.SWR_modulation;
unfiltered_SWR_modulation = metrics.unfiltered_SWR_modulation;
SWR_modulation_real_duration = metrics.SWR_modulation_real_duration;
swr_duration = metrics.swr_duration;
SWR_all_percentage_reactivation = metrics.SWR_all_percentage_reactivation;
SWR_all_binary = metrics.SWR_all_binary;
peaks_amp_sd = metrics.peaks_amp_sd;
PSP = metrics.PSP;
SFI = metrics.SFI;
SSMI = metrics.SSMI;
longPSP = metrics.longPSP;
SpPR = metrics.SpPR;
SRR = metrics.SRR;

[swr_data_dir, ~, ~] = fileparts(swr_data_file);
if exist(swr_data_dir, 'dir') ~= 7
    mkdir(swr_data_dir);
end
save(swr_data_file, 'swr_duration_mean', 'SWR_spiking_activity', ...
    'SWR_activation_durations', 'SWR_number', 'non_zero_mean_first_time_of_spikes', ...
    'non_zero_median_spike_time', 'non_zero_mean_spike_time', 'time_of_spikes', ...
    'SWRspikes', 'SWRrealspikes', 'pre_SWR_real_spikes', 'pre_SWRspikes', ...
    'SWR_sums', 'SWR_real_sums', 'pre_SWR_real_sums', 'pre_SWR_sums', ...
    'SWR_modulation', 'unfiltered_SWR_modulation', 'SWR_modulation_real_duration', ...
    'swr_duration', 'SWR_all_percentage_reactivation', 'SWR_all_binary', ...
    'peaks_amp_sd', 'PSP', 'SFI', 'SSMI', 'longPSP', 'SpPR', 'SRR');
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

function tf = window_overlaps_events(window_start, window_stop, event_starts, event_stops)

tf = false;
for event_idx = 1:numel(event_starts)
    if (window_start >= event_starts(event_idx) && window_start <= event_stops(event_idx)) || ...
            (window_stop >= event_starts(event_idx) && window_stop <= event_stops(event_idx)) || ...
            (window_start <= event_starts(event_idx) && window_stop >= event_stops(event_idx))
        tf = true;
        return
    end
end
end

function duration_sum = sum_swr_duration_field(slow_cHFOs, swr_duration)

if isfield(slow_cHFOs, 'duration')
    duration_sum = sum([slow_cHFOs.duration]);
else
    duration_sum = sum(swr_duration);
end
end

function values = column_nanmean(matrix_values)

if isempty(matrix_values)
    values = [];
    return
end

valid_counts = sum(isfinite(matrix_values), 1);
matrix_values(~isfinite(matrix_values)) = 0;
values = sum(matrix_values, 1) ./ valid_counts;
values(valid_counts == 0) = NaN;
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

function metrics = empty_metrics()

metrics = struct( ...
    'status', '', ...
    'sleepLabel', '', ...
    'swrDataFile', '', ...
    'SWR_number', 0, ...
    'swr_duration_mean', NaN, ...
    'SWR_spiking_activity', {{}}, ...
    'SWR_activation_durations', {{}}, ...
    'non_zero_mean_first_time_of_spikes', [], ...
    'non_zero_median_spike_time', [], ...
    'non_zero_mean_spike_time', [], ...
    'time_of_spikes', {{}}, ...
    'SWRspikes', [], ...
    'SWRrealspikes', [], ...
    'pre_SWR_real_spikes', [], ...
    'pre_SWRspikes', [], ...
    'SWR_sums', [], ...
    'SWR_real_sums', [], ...
    'pre_SWR_real_sums', [], ...
    'pre_SWR_sums', [], ...
    'SWR_modulation', [], ...
    'unfiltered_SWR_modulation', [], ...
    'SWR_modulation_real_duration', [], ...
    'swr_duration', [], ...
    'SWR_all_percentage_reactivation', [], ...
    'SWR_all_binary', [], ...
    'peaks_amp_sd', [], ...
    'PSP', [], ...
    'SFI', [], ...
    'SSMI', [], ...
    'longPSP', [], ...
    'SpPR', [], ...
    'SRR', [] ...
    );
end

function summary = empty_swr_summary()

summary = struct( ...
    'session_index', NaN, ...
    'status', '', ...
    'cell_count', 0, ...
    'processed_sleeps', {{}}, ...
    'sleep_statuses', {{}}, ...
    'swr_counts', [] ...
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
