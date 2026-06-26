function [slow_cHFOs, fast_cHFOs, cHFOs, lfp_bp] = detect_swr_events_blanco(data, filter_numerator, detection_threshold, velocities)
% Detect candidate SWR/HFO events with the Blanco-style RMS threshold.

lfp = data.samp;
sample_frequency = data.sampFreq(1);

detection_threshold = fill_detection_defaults(detection_threshold);

lfp_bp = filtfilt(filter_numerator, 1, lfp);

remainder = mod(numel(lfp_bp), 6) - 1;
if remainder >= 0
    lfp_bp(end - remainder:end) = [];
end

RMS = zeros(size(lfp_bp));
for sample_idx = 1:numel(lfp_bp) - 5
    window = lfp_bp(sample_idx:sample_idx + 5);
    RMS(sample_idx) = sqrt(mean(window .^ 2));
end

RMS_sd = std(RMS);
RMS_mean = mean(RMS);
threshold = RMS_mean + detection_threshold.numSD * RMS_sd;
logical_RMS = RMS > threshold;

[start_idx, stop_idx] = detectHFOs.start_stop(logical_RMS);
if isempty(start_idx)
    [slow_cHFOs, fast_cHFOs, cHFOs] = empty_event_outputs(lfp_bp);
    return
end

event_duration_samples = stop_idx - start_idx;
too_short = event_duration_samples < detection_threshold.minDur * sample_frequency;
start_idx(too_short) = [];
stop_idx(too_short) = [];

event_duration_samples = stop_idx - start_idx;
too_long = event_duration_samples > detection_threshold.maxDur * sample_frequency;
start_idx(too_long) = [];
stop_idx(too_long) = [];

if isempty(start_idx)
    [slow_cHFOs, fast_cHFOs, cHFOs] = empty_event_outputs(lfp_bp);
    return
end

interval_threshold = round(detection_threshold.mergeThreshold * sample_frequency);
merge_next = zeros(numel(start_idx) - 1, 1);
for event_idx = 1:numel(start_idx) - 1
    merge_next(event_idx) = (start_idx(event_idx + 1) - stop_idx(event_idx)) < interval_threshold;
end
merge_next = logical(merge_next);
start_idx([false; merge_next]) = [];
stop_idx([merge_next; false]) = [];

if isempty(start_idx)
    [slow_cHFOs, fast_cHFOs, cHFOs] = empty_event_outputs(lfp_bp);
    return
end

rectified_lfp_bp = abs(lfp_bp);
peak_threshold = detection_threshold.peakSD * std(rectified_lfp_bp);
discard_event = false(size(start_idx));
peak_amp_sd = zeros(size(start_idx));

for event_idx = 1:numel(start_idx)
    event = rectified_lfp_bp(start_idx(event_idx):stop_idx(event_idx));
    peaks = findpeaks(event);
    if isempty(peaks)
        peak_amp_sd(event_idx) = NaN;
    else
        peak_amp_sd(event_idx) = (max(peaks) - RMS_mean) / RMS_sd;
    end
    peak_locations = findpeaksmine(event);
    peak_locations = peak_locations.loc;
    discard_event(event_idx) = sum(event(peak_locations) > peak_threshold) < detection_threshold.numHighPwrCycles * 2;
end

start_idx(discard_event) = [];
stop_idx(discard_event) = [];
peak_amp_sd(discard_event) = [];

if isempty(start_idx)
    [slow_cHFOs, fast_cHFOs, cHFOs] = empty_event_outputs(lfp_bp);
    return
end

cHFOs = build_event_structs(data, start_idx, stop_idx, peak_amp_sd);
slow_mask = find_slow_event_mask(data, cHFOs, velocities);
slow_cHFOs = cHFOs(slow_mask);
fast_cHFOs = cHFOs(~slow_mask);
end

function detection_threshold = fill_detection_defaults(detection_threshold)

if nargin < 1 || isempty(detection_threshold)
    detection_threshold = struct();
end
if ~isfield(detection_threshold, 'numSD')
    detection_threshold.numSD = 3.5;
end
if ~isfield(detection_threshold, 'minDur')
    detection_threshold.minDur = 0.006;
end
if ~isfield(detection_threshold, 'maxDur')
    detection_threshold.maxDur = 0.5;
end
if ~isfield(detection_threshold, 'numCycles')
    detection_threshold.numCycles = 6;
end
if ~isfield(detection_threshold, 'numHighPwrCycles')
    detection_threshold.numHighPwrCycles = 5;
end
if ~isfield(detection_threshold, 'mergeThreshold')
    detection_threshold.mergeThreshold = 0.01;
end
if ~isfield(detection_threshold, 'peakSD')
    detection_threshold.peakSD = 3;
end
end

function cHFOs = build_event_structs(data, start_idx, stop_idx, peak_amp_sd)

cHFOs = repmat(empty_event_struct(), 1, numel(start_idx));
for event_idx = 1:numel(start_idx)
    duration = data.ts(stop_idx(event_idx)) - data.ts(start_idx(event_idx));
    cHFOs(event_idx) = struct( ...
        'start_ind', start_idx(event_idx), ...
        'stop_ind', stop_idx(event_idx), ...
        'start_ts', data.ts(start_idx(event_idx)), ...
        'stop_ts', data.ts(stop_idx(event_idx)), ...
        'duration', duration, ...
        'peaks_sd', peak_amp_sd(event_idx), ...
        'ttNo', 0, ...
        'cluster', 0, ...
        'velocity', NaN, ...
        'x', NaN, ...
        'y', NaN, ...
        'princ_freq', NaN, ...
        'princ_freq_amp', NaN, ...
        'mean_amp', NaN, ...
        'spikes_ts', cell(1, 1));
end
end

function slow_mask = find_slow_event_mask(data, cHFOs, velocities)

slow_mask = false(1, numel(cHFOs));
slow_velocities = velocities < 2;
for event_idx = 1:numel(cHFOs)
    swr_mid_point = cHFOs(event_idx).start_ts + ...
        (cHFOs(event_idx).stop_ts - cHFOs(event_idx).start_ts) / 2;

    before = swr_mid_point - 1.5;
    after = swr_mid_point + 1.5;
    [~, closest_idx_before] = min(abs(data.ts - before));
    [~, closest_idx_after] = min(abs(data.ts - after));
    start_window_idx = min(closest_idx_before, closest_idx_after);
    stop_window_idx = max(closest_idx_before, closest_idx_after);

    slow_mask(event_idx) = all(slow_velocities(start_window_idx:stop_window_idx) == 1);
end
end

function [slow_cHFOs, fast_cHFOs, cHFOs] = empty_event_outputs(~)

cHFOs = repmat(empty_event_struct(), 0, 0);
slow_cHFOs = cHFOs;
fast_cHFOs = cHFOs;
end

function event = empty_event_struct()

event = struct( ...
    'start_ind', [], ...
    'stop_ind', [], ...
    'start_ts', [], ...
    'stop_ts', [], ...
    'duration', [], ...
    'peaks_sd', [], ...
    'ttNo', [], ...
    'cluster', [], ...
    'velocity', [], ...
    'x', [], ...
    'y', [], ...
    'princ_freq', [], ...
    'princ_freq_amp', [], ...
    'mean_amp', [], ...
    'spikes_ts', []);
end
