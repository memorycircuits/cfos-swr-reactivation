function SALT_opto_analysis(custom_settings)
% SALT_opto_analysis
% Compute SALT-style optotagging metrics and classify tagged cells.
%
% Direct and indirect tagged cells are separated with the settled fixed
% latency rule: latencies below 7 ms are direct; latencies at or above
% 7 ms are indirect. This file writes the optotagged field, saves the
% classification cache, and generates figures tied directly to the SALT
% procedure. CSV exports are handled separately by CSV_file_export_first.

close all
clc
rng(1)

global SALT_opto_analysis_custom_settings

if nargin < 1 || isempty(custom_settings)
    if ~isempty(SALT_opto_analysis_custom_settings)
        custom_settings = SALT_opto_analysis_custom_settings;
    else
        custom_settings = struct();
    end
end

preferredOutDir = get_override_value(custom_settings, 'preferredOutDir', fullfile(pwd, 'OptoMetricComparison_SALT'));
analyze_dir = char(string(get_override_value(custom_settings, 'analyze_dir', 'opto1')));
default_salt_cutoff = get_override_value(custom_settings, 'default_salt_cutoff', 0.01);
min_reliability_for_tagged = get_override_value(custom_settings, 'min_reliability_for_tagged', 0.10);
reuse_existing_salt = get_override_value(custom_settings, 'reuse_existing_salt', true);
principal_codes = get_override_value(custom_settings, 'principal_codes', 1);
interneuron_codes = get_override_value(custom_settings, 'interneuron_codes', 2);
response_window_ms = get_override_value(custom_settings, 'response_window_ms', 15);
delay_ms = get_override_value(custom_settings, 'delay_ms', 0);
guard_ms = get_override_value(custom_settings, 'guard_ms', 5);
latency_bin_ms = get_override_value(custom_settings, 'latency_bin_ms', 1);
baseline_bootstrap_sets = get_override_value(custom_settings, 'baseline_bootstrap_sets', 100);
null_pair_draws = get_override_value(custom_settings, 'null_pair_draws', 500);
min_baseline_windows = get_override_value(custom_settings, 'min_baseline_windows', 20);
fixed_latency_threshold_ms = 7;
writeback_field_name = char(string(get_override_value(custom_settings, 'writeback_field_name', 'optotagged')));

stageB_settings = ensure_stage_b_settings_defaults(struct('fixed_latency_threshold_ms', fixed_latency_threshold_ms));

response_window_s = response_window_ms / 1000;
delay_s = delay_ms / 1000;
guard_s = guard_ms / 1000;
latency_bin_s = latency_bin_ms / 1000;

outDir = resolve_output_dir(preferredOutDir, stageB_settings);
allCellsPath = char(string(get_override_value(custom_settings, 'allCellsPath', '')));
if isempty(allCellsPath)
    allCellsPath = resolve_all_cells_path();
end
sessionInfoPath = char(string(get_override_value(custom_settings, 'sessionInfoPath', '')));
if isempty(sessionInfoPath)
    sessionInfoPath = resolve_session_info_path();
end
resultsCachePath = fullfile(outDir, 'SALT_TwoStep_ClassificationResults.mat');

fprintf('\n=== SALT Optotagging Analysis ===\n');
fprintf('Stage B rule: %s\n', describe_stage_b_export_variant(stageB_settings));
fprintf('Results directory: %s\n', outDir);

if ~exist(outDir, 'dir')
    mkdir(outDir);
end

add_analysis_paths(custom_settings);

loadedAllCells = load(allCellsPath, 'All_Cells_combined');
All_Cells_combined = loadedAllCells.All_Cells_combined;
availableVariables = who('-file', allCellsPath);
if any(strcmp(availableVariables, 'last_updated'))
    loadedLastUpdated = load(allCellsPath, 'last_updated');
    last_updated = loadedLastUpdated.last_updated;
    for iSession = 1:numel(last_updated)
        sessionFields = fieldnames(last_updated(iSession));
        for iField = 1:numel(sessionFields)
            fieldName = sessionFields{iField};
            fieldValue = last_updated(iSession).(fieldName);

            if isstruct(fieldValue)
                if ~isfield(fieldValue, 'update_timestamp')
                    fieldValue.update_timestamp = '';
                end
                if ~isfield(fieldValue, 'update_script')
                    fieldValue.update_script = 'legacy_unknown';
                end
            else
                if ischar(fieldValue)
                    timestampValue = fieldValue;
                elseif isstring(fieldValue) && isscalar(fieldValue)
                    timestampValue = char(fieldValue);
                else
                    timestampValue = '';
                end
                fieldValue = struct( ...
                    'update_timestamp', timestampValue, ...
                    'update_script', 'legacy_unknown');
            end

            last_updated(iSession).(fieldName) = fieldValue;
        end
    end
else
    last_updated = struct;
end

load(sessionInfoPath, 'sessInfo');

saltPBySession = cell(numel(All_Cells_combined), 1);
saltZBySession = cell(numel(All_Cells_combined), 1);
saltStatsBySession = cell(numel(All_Cells_combined), 1);

validSessions = find(arrayfun(@(s) max_num_cells_in_session(s) > 0, All_Cells_combined));
validSessions = validSessions(validSessions <= numel(sessInfo));

processedSessions = 0;
skippedSessions = 0;

for iSession = 1:numel(validSessions)
    s = validSessions(iSession);
    fprintf('\nSession %d/%d (i=%d)\n', iSession, numel(validSessions), s);

    if ~isfield(sessInfo(s), 'mainDir') || ~isfield(sessInfo(s), 'tList')
        fprintf('  Missing session metadata. Skipping.\n');
        skippedSessions = skippedSessions + 1;
        continue
    end

    mainDir = normalize_text_field(sessInfo(s).mainDir);
    tListName = normalize_text_field(sessInfo(s).tList);

    if isempty(mainDir) || isempty(tListName)
        fprintf('  Missing or malformed mainDir/tList. Skipping.\n');
        skippedSessions = skippedSessions + 1;
        continue
    end

    tListPath = fullfile(mainDir, tListName);
    optoDir = fullfile(mainDir, analyze_dir);

    if ~exist(tListPath, 'file')
        fprintf('  Missing tList file. Skipping.\n');
        skippedSessions = skippedSessions + 1;
        continue
    end

    if ~exist(optoDir, 'dir')
        fprintf('  Missing %s folder. Skipping.\n', analyze_dir);
        skippedSessions = skippedSessions + 1;
        continue
    end

    tt_files = read_tt_list(tListPath);
    nCells = numel(tt_files);
    if nCells == 0
        fprintf('  No cells found in tList. Skipping.\n');
        skippedSessions = skippedSessions + 1;
        continue
    end

    saltMetricPath = fullfile(optoDir, sprintf('SALTMetric_%dmsDelay.mat', delay_ms));
    if reuse_existing_salt && exist(saltMetricPath, 'file')
        [cacheOk, cachedMetrics, cacheMessage] = try_load_cached_salt_metrics( ...
            saltMetricPath, nCells, response_window_ms, delay_ms, guard_ms, ...
            latency_bin_ms, baseline_bootstrap_sets, null_pair_draws, ...
            min_baseline_windows, analyze_dir);

        if cacheOk
            saltPBySession{s} = cachedMetrics.salt_p_value(:);
            saltZBySession{s} = cachedMetrics.salt_z_value(:);
            saltStatsBySession{s} = cachedMetrics.salt_statistic(:);

            fprintf('  %s\n', cacheMessage);
            fprintf('  Cells: %d\n', nCells);
            fprintf('  SALT-style significant (p < %.3f): %d\n', default_salt_cutoff, ...
                nnz(cachedMetrics.salt_p_value(:) < default_salt_cutoff));

            processedSessions = processedSessions + 1;
            continue
        elseif ~isempty(cacheMessage)
            fprintf('  %s\n', cacheMessage);
        end
    end

    try
        [TTLData.ON_ts, TTLData.OFF_ts, TTLData.lengthONStim] = GetLaserTTL_hj(optoDir); %#ok<NASGU>
    catch ME
        fprintf('  Could not load TTLs: %s\n', ME.message);
        skippedSessions = skippedSessions + 1;
        continue
    end

    if isempty(TTLData.ON_ts) || isempty(TTLData.OFF_ts)
        fprintf('  No TTL pulses found. Skipping.\n');
        skippedSessions = skippedSessions + 1;
        continue
    end

    try
        spikeData = readSpikeDataOnly(optoDir, tt_files);
        spikeTimes = fixSpikes(spikeData);
    catch ME
        fprintf('  Could not load spike data: %s\n', ME.message);
        skippedSessions = skippedSessions + 1;
        continue
    end

    ttlOn = TTLData.ON_ts(:) + delay_s;
    ttlOff = TTLData.OFF_ts(:);
    nPulses = numel(ttlOn);
    baselineWindows = build_baseline_windows(ttlOn, ttlOff, response_window_s, guard_s);
    nBaseline = size(baselineWindows, 1);

    if nPulses == 0 || nBaseline < min_baseline_windows
        fprintf('  Not enough pulses or baseline windows for SALT. Skipping.\n');
        skippedSessions = skippedSessions + 1;
        continue
    end

    nLatencyBins = max(1, ceil(response_window_s / latency_bin_s));
    nCategories = nLatencyBins + 1;

    baselineSampleIdx = generate_baseline_sample_idx(nBaseline, nPulses, baseline_bootstrap_sets);
    nullPairIdx = generate_null_pair_idx(baseline_bootstrap_sets, null_pair_draws);

    pVals = nan(nCells, 1);
    zVals = nan(nCells, 1);
    statVals = nan(nCells, 1);

    for iCell = 1:nCells
        if iCell > numel(spikeTimes) || isempty(spikeTimes{iCell})
            continue
        end

        cellSpikes = spikeTimes{iCell}(:);
        stimCats = first_spike_categories(cellSpikes, ttlOn, response_window_s, latency_bin_s, nLatencyBins);
        baseCats = first_spike_categories(cellSpikes, baselineWindows(:,1), response_window_s, latency_bin_s, nLatencyBins);

        [pVals(iCell), zVals(iCell), statVals(iCell)] = compute_salt_style_metric( ...
            stimCats, baseCats, baselineSampleIdx, nullPairIdx, nCategories);
    end

    saltPBySession{s} = pVals;
    saltZBySession{s} = zVals;
    saltStatsBySession{s} = statVals;

    sessionMetrics = struct();
    sessionMetrics.parameters = struct( ...
        'response_window_ms', response_window_ms, ...
        'delay_ms', delay_ms, ...
        'guard_ms', guard_ms, ...
        'latency_bin_ms', latency_bin_ms, ...
        'baseline_bootstrap_sets', baseline_bootstrap_sets, ...
        'null_pair_draws', null_pair_draws, ...
        'min_baseline_windows', min_baseline_windows, ...
        'analyze_dir', analyze_dir);
    sessionMetrics.salt_p_value = pVals;
    sessionMetrics.salt_z_value = zVals;
    sessionMetrics.salt_statistic = statVals;
    sessionMetrics.tt_files = tt_files;
    sessionMetrics.ttl_on = ttlOn;
    sessionMetrics.ttl_off = ttlOff;
    sessionMetrics.baseline_windows = baselineWindows;

    save(saltMetricPath, 'sessionMetrics');

    fprintf('  Cells: %d\n', nCells);
    fprintf('  Pulses: %d\n', nPulses);
    fprintf('  Baseline windows: %d\n', nBaseline);
    fprintf('  SALT-style significant (p < %.3f): %d\n', default_salt_cutoff, nnz(pVals < default_salt_cutoff));

    processedSessions = processedSessions + 1;
end

fprintf('\nProcessed sessions: %d\n', processedSessions);
fprintf('Skipped sessions: %d\n', skippedSessions);

T = flatten_cells_for_salt(All_Cells_combined, saltPBySession, principal_codes, interneuron_codes);
T.Latency_ms(T.Reliability <= 0) = NaN;
T.Jitter_ms(T.Reliability <= 0) = NaN;
T.Latency_ms(T.Latency_ms < 0) = NaN;
T.Jitter_ms(T.Jitter_ms < 0) = NaN;
[T, SaltTwoStep] = run_salt_two_step_classification( ...
    T, default_salt_cutoff, min_reliability_for_tagged, stageB_settings);

[All_Cells_combined, last_updated, optoCounts] = write_optotagged_field_to_all_cells( ...
    All_Cells_combined, last_updated, T, writeback_field_name, 'SALT_opto_analysis.m');
save(allCellsPath, 'All_Cells_combined', 'last_updated', '-v7.3');

ExportSettings = struct( ...
    'default_salt_cutoff', default_salt_cutoff, ...
    'min_reliability_for_tagged', min_reliability_for_tagged, ...
    'stageB_settings', stageB_settings, ...
    'analyze_dir', analyze_dir, ...
    'delay_ms', delay_ms, ...
    'response_window_ms', response_window_ms, ...
    'guard_ms', guard_ms, ...
    'writeback_field_name', writeback_field_name);
save(resultsCachePath, 'T', 'SaltTwoStep', 'ExportSettings', '-v7.3');
run_salt_procedure_figures(T, SaltTwoStep, ExportSettings, outDir);

fprintf('\nSaved updated All_Cells_combined\n');
fprintf('Field written: %s\n', writeback_field_name);
fprintf('Coding: direct = 1, not responding = 0, indirect = -1\n');
fprintf('Stage A rule: SALT p < %.4g and reliability >= %.3f\n', default_salt_cutoff, min_reliability_for_tagged);
fprintf('Stage B rule: direct latency < %.3g ms; indirect latency >= %.3g ms\n', fixed_latency_threshold_ms, fixed_latency_threshold_ms);
fprintf('Final counts written:\n');
fprintf('  Direct: %d\n', optoCounts.Direct);
fprintf('  Not responding: %d\n', optoCounts.NotResponding);
fprintf('  Indirect: %d\n', optoCounts.Indirect);
fprintf('\nSaved SALT classification cache to:\n%s\n', resultsCachePath);

end

function nCells = max_num_cells_in_session(S)
    fieldsToCheck = {'SD_0_ms_delay', 'spike_probability_per_TTL', ...
        'cells_mode_TTL_latency', 'cells_std_TTL_latency', ...
        'GMM_based_classification_days'};

    nCells = 0;
    for i = 1:numel(fieldsToCheck)
        if isfield(S, fieldsToCheck{i})
            x = S.(fieldsToCheck{i});
            if ~isempty(x)
                nCells = max(nCells, numel(x));
            end
        end
    end
end

function [All_Cells_combined, last_updated, counts] = write_optotagged_field_to_all_cells(All_Cells_combined, last_updated, T, field_name, update_script_name)

    counts = struct('Direct', 0, 'NotResponding', 0, 'Indirect', 0);
    update_timestamp = datestr(now, 'yyyy-mm-ddTHH:MM:SS');

    for s = 1:numel(All_Cells_combined)
        nCells = max_num_cells_in_session(All_Cells_combined(s));
        if nCells == 0
            All_Cells_combined(s).(field_name) = [];
            last_updated = mark_field_updated(last_updated, s, field_name, update_timestamp, update_script_name);
            continue
        end

        optotagged = zeros(nCells, 1);
        sessionMask = T.SessionIndex == s;

        if any(sessionMask)
            cellIdx = T.Cell(sessionMask);
            finalLabel = T.FinalLabel(sessionMask);
            validIdx = isfinite(cellIdx) & cellIdx >= 1 & cellIdx <= nCells;

            cellIdx = cellIdx(validIdx);
            finalLabel = finalLabel(validIdx);

            directMask = finalLabel == "direct";
            indirectMask = finalLabel == "indirect";

            optotagged(cellIdx(directMask)) = 1;
            optotagged(cellIdx(indirectMask)) = -1;
        end

        All_Cells_combined(s).(field_name) = optotagged;
        last_updated = mark_field_updated(last_updated, s, field_name, update_timestamp, update_script_name);
        counts.Direct = counts.Direct + nnz(optotagged == 1);
        counts.NotResponding = counts.NotResponding + nnz(optotagged == 0);
        counts.Indirect = counts.Indirect + nnz(optotagged == -1);
    end
end

function idxMat = generate_baseline_sample_idx(nBaseline, nPulses, nSets)

    idxMat = zeros(nPulses, nSets);
    sampleWithReplacement = nBaseline < nPulses;

    for i = 1:nSets
        if sampleWithReplacement
            idxMat(:, i) = randi(nBaseline, nPulses, 1);
        else
            idxMat(:, i) = randperm(nBaseline, nPulses)';
        end
    end
end

function pairIdx = generate_null_pair_idx(nSets, nPairs)

    pairIdx = zeros(nPairs, 2);
    for i = 1:nPairs
        a = randi(nSets);
        b = randi(nSets);
        while b == a
            b = randi(nSets);
        end
        pairIdx(i,:) = [a b];
    end
end

function cats = first_spike_categories(spikeTimes, windowStarts, response_window_s, latency_bin_s, nLatencyBins)

    nWindows = numel(windowStarts);
    cats = repmat(nLatencyBins + 1, nWindows, 1);

    for i = 1:nWindows
        wStart = windowStarts(i);
        wEnd = wStart + response_window_s;
        idx = find(spikeTimes >= wStart & spikeTimes < wEnd, 1, 'first');
        if isempty(idx)
            continue
        end

        latency = spikeTimes(idx) - wStart;
        bin = floor(latency / latency_bin_s) + 1;
        bin = max(1, min(bin, nLatencyBins));
        cats(i) = bin;
    end
end

function [pVal, zVal, testStat] = compute_salt_style_metric(stimCats, baseCats, baselineSampleIdx, nullPairIdx, nCategories)

    if isempty(stimCats) || isempty(baseCats)
        pVal = NaN;
        zVal = NaN;
        testStat = NaN;
        return
    end

    histStim = category_histogram(stimCats, nCategories);
    nSets = size(baselineSampleIdx, 2);
    baselineHists = zeros(nSets, nCategories);

    for i = 1:nSets
        baselineHists(i,:) = category_histogram(baseCats(baselineSampleIdx(:, i)), nCategories);
    end

    testDistances = zeros(nSets, 1);
    for i = 1:nSets
        testDistances(i) = jensen_shannon_divergence(histStim, baselineHists(i,:));
    end

    nullDistances = zeros(size(nullPairIdx, 1), 1);
    for i = 1:size(nullPairIdx, 1)
        nullDistances(i) = jensen_shannon_divergence( ...
            baselineHists(nullPairIdx(i,1),:), baselineHists(nullPairIdx(i,2),:));
    end

    testStat = median(testDistances);
    pVal = (1 + nnz(nullDistances >= testStat)) / (numel(nullDistances) + 1);

    nullMean = mean(nullDistances, 'omitnan');
    nullStd = std(nullDistances, 0, 'omitnan');
    if isnan(nullStd) || nullStd == 0
        zVal = NaN;
    else
        zVal = (testStat - nullMean) / nullStd;
    end
end

function T = flatten_cells_for_salt(All_Cells_combined, saltPBySession, principal_codes, interneuron_codes)

    SessionIndex = [];
    INumber = [];
    Cell = [];
    CellClass = strings(0,1);
    AnimalID = strings(0,1);
    RoomID = [];
    RoomName = strings(0,1);
    RoomGroup = strings(0,1);
    BaselineSD = [];
    Reliability = [];
    Latency_ms = [];
    Jitter_ms = [];
    SALTPValue = [];
    S1_SSMI = [];
    S1_SFI = [];
    S1_PSP = [];
    S1_SpPR = [];
    S2_SSMI = [];
    S2_SFI = [];
    S2_PSP = [];
    S2_SpPR = [];
    BurstIndex = [];
    SpeedScore = [];
    PeakRate = [];
    TemporalPeakRate = [];
    AverageRate = [];
    AverageRate_OF1 = [];
    AverageRate_OF2 = [];
    AverageRate_OF3 = [];
    PlaceFieldSize = [];
    SpatialInfo_A = [];
    SpatialPeakRate = [];
    SpatialPeakRate_OF1 = [];
    SpatialPeakRate_OF2 = [];
    SpatialPeakRate_OF3 = [];
    PlaceFieldSize_OF1 = [];
    PlaceFieldSize_OF2 = [];
    PlaceFieldSize_OF3 = [];
    PlaceFieldNumber = [];
    PlaceFieldNumber_OF1 = [];
    PlaceFieldNumber_OF2 = [];
    PlaceFieldNumber_OF3 = [];
    ZCoherence_OF1 = [];
    ZCoherence_OF2 = [];
    ZCoherence_OF3 = [];
    Sparseness_OF1 = [];
    Sparseness_OF2 = [];
    Sparseness_OF3 = [];
    Selectivity_OF1 = [];
    Selectivity_OF2 = [];
    Selectivity_OF3 = [];
    SpatialInfo_OF1 = [];
    SpatialInfo_OF2 = [];
    SpatialInfo_OF3 = [];
    SpatialInfoRate_OF1 = [];
    SpatialInfoRate_OF2 = [];
    SpatialInfoRate_OF3 = [];
    Remapping = [];
    Stability = [];
    ThetaModulationScore = [];
    ThetaModulationTanaka = [];
    IntrinsicFrequency = [];
    ThetaPhasePreferred = [];
    ThetaPhaseR = [];
    ThetaPhaseP = [];
    ThetaPhaseZ = [];
    ThetaPhaseNSpikes = [];
    AverageRate_OFMean = [];
    PeakRate_OFMax = [];
    SpatialPeakRate_OFMax = [];
    TemporalPeakRate_OF1 = [];
    TemporalPeakRate_OF2 = [];
    TemporalPeakRate_OF3 = [];
    BurstIndex_OFMean = [];
    ThetaModulationScore_OFMean = [];
    ThetaModulationTanaka_OFMean = [];
    IntrinsicFrequency_OFMean = [];

    for s = 1:numel(All_Cells_combined)
        bm = get_field_or_empty(All_Cells_combined(s), 'SD_0_ms_delay');
        rel = get_field_or_empty(All_Cells_combined(s), 'spike_probability_per_TTL');
        lat = get_field_or_empty(All_Cells_combined(s), 'cells_mode_TTL_latency');
        jit = get_field_or_empty(All_Cells_combined(s), 'cells_std_TTL_latency');
        cls = get_field_or_empty(All_Cells_combined(s), 'GMM_based_classification_days');
        animal = get_field_or_empty(All_Cells_combined(s), 'animal');
        iNumber = get_field_or_empty(All_Cells_combined(s), 'i_number');
        s1_ssmi = get_field_or_empty(All_Cells_combined(s), 'S1_SSMI');
        s1_sfi = get_field_or_empty(All_Cells_combined(s), 'S1_SFI');
        s1_psp = get_field_or_empty(All_Cells_combined(s), 'S1_PSP');
        s1_sppr = get_field_or_empty(All_Cells_combined(s), 'S1_SpPR');
        s2_ssmi = get_field_or_empty(All_Cells_combined(s), 'S2_SSMI');
        s2_sfi = get_field_or_empty(All_Cells_combined(s), 'S2_SFI');
        s2_psp = get_field_or_empty(All_Cells_combined(s), 'S2_PSP');
        s2_sppr = get_field_or_empty(All_Cells_combined(s), 'S2_SpPR');
        room = get_field_or_empty(All_Cells_combined(s), 'room_ID');
        burst = get_field_or_empty(All_Cells_combined(s), 'of_avg_burst_indices');
        speedScore = get_field_or_empty(All_Cells_combined(s), 'of_avg_speedScore');
        peakRate = get_first_available_field(All_Cells_combined(s), {'peak_rate', 'of_avg_peak_fr_rate'});
        avgRate = get_field_or_empty(All_Cells_combined(s), 'classific_firingRate');
        avgRateOF = get_first_available_field(All_Cells_combined(s), {'of_avg_fir_rate', 'classific_firingRate'});
        avgRate1 = get_first_available_field(All_Cells_combined(s), {'averRate_of1', 'AverageRate_OF1', 'meanRate_of1', 'of1_fir_rate', 'fir_rate_of1'});
        avgRate2 = get_first_available_field(All_Cells_combined(s), {'averRate_of2', 'AverageRate_OF2', 'meanRate_of2', 'of2_fir_rate', 'fir_rate_of2'});
        avgRate3 = get_first_available_field(All_Cells_combined(s), {'averRate_of3', 'AverageRate_OF3', 'meanRate_of3', 'of3_fir_rate', 'fir_rate_of3'});
        peakRate1 = get_first_available_field(All_Cells_combined(s), {'spatialPeakRate_of1', 'SpatialPeakRate_OF1', 'peakRate_of1'});
        peakRate2 = get_first_available_field(All_Cells_combined(s), {'spatialPeakRate_of2', 'SpatialPeakRate_OF2', 'peakRate_of2'});
        peakRate3 = get_first_available_field(All_Cells_combined(s), {'spatialPeakRate_of3', 'SpatialPeakRate_OF3', 'peakRate_of3'});
        temporalPeakRate1 = get_temporal_peak_rate_of(All_Cells_combined(s), 1);
        temporalPeakRate2 = get_temporal_peak_rate_of(All_Cells_combined(s), 2);
        temporalPeakRate3 = get_temporal_peak_rate_of(All_Cells_combined(s), 3);
        burst1 = get_of_metric_or_empty(All_Cells_combined(s), ...
            {'burst_index_of1', 'burstIndex_of1', 'of1_burst_index', 'of1_burstIndex', 'burst_indices_of1'}, ...
            {'burst_index_sess', 'burstIndex_sess', 'burst_indices_sess'}, 1);
        burst2 = get_of_metric_or_empty(All_Cells_combined(s), ...
            {'burst_index_of2', 'burstIndex_of2', 'of2_burst_index', 'of2_burstIndex', 'burst_indices_of2'}, ...
            {'burst_index_sess', 'burstIndex_sess', 'burst_indices_sess'}, 2);
        burst3 = get_of_metric_or_empty(All_Cells_combined(s), ...
            {'burst_index_of3', 'burstIndex_of3', 'of3_burst_index', 'of3_burstIndex', 'burst_indices_of3'}, ...
            {'burst_index_sess', 'burstIndex_sess', 'burst_indices_sess'}, 3);
        pf1 = get_first_available_field(All_Cells_combined(s), {'PF_sumSize_of1', 'context_A_field_size'});
        pf2 = get_first_available_field(All_Cells_combined(s), {'PF_sumSize_of2', 'context_B_field_size'});
        pf3 = get_first_available_field(All_Cells_combined(s), {'PF_sumSize_of3', 'context_A_revisit_field_size'});
        pfnum1 = get_first_available_field(All_Cells_combined(s), {'PF_fieldNumbers_of1', 'of1_place_field_numbers'});
        pfnum2 = get_first_available_field(All_Cells_combined(s), {'PF_fieldNumbers_of2', 'of2_place_field_numbers'});
        pfnum3 = get_first_available_field(All_Cells_combined(s), {'PF_fieldNumbers_of3', 'of3_place_field_numbers'});
        zcoh1 = get_first_available_field(All_Cells_combined(s), {'PF_zCoherence_of1'});
        zcoh2 = get_first_available_field(All_Cells_combined(s), {'PF_zCoherence_of2'});
        zcoh3 = get_first_available_field(All_Cells_combined(s), {'PF_zCoherence_of3'});
        sparse1 = get_first_available_field(All_Cells_combined(s), {'PF_sparseness_of1'});
        sparse2 = get_first_available_field(All_Cells_combined(s), {'PF_sparseness_of2'});
        sparse3 = get_first_available_field(All_Cells_combined(s), {'PF_sparseness_of3'});
        select1 = get_first_available_field(All_Cells_combined(s), {'PF_selectivity_of1'});
        select2 = get_first_available_field(All_Cells_combined(s), {'PF_selectivity_of2'});
        select3 = get_first_available_field(All_Cells_combined(s), {'PF_selectivity_of3'});
        spatialInfo_A = get_field_or_empty(All_Cells_combined(s), 'context_A_spatial_info');
        spatialInfo_B = get_field_or_empty(All_Cells_combined(s), 'context_B_spatial_info');
        spatialInfo_A_revisit = get_field_or_empty(All_Cells_combined(s), 'context_A_revisit_spatial_info');
        spatialInfoRate_A = get_field_or_empty(All_Cells_combined(s), 'context_A_spatial_info_rate');
        spatialInfoRate_B = get_field_or_empty(All_Cells_combined(s), 'context_B_spatial_info_rate');
        spatialInfoRate_A_revisit = get_field_or_empty(All_Cells_combined(s), 'context_A_revisit_spatial_info_rate');
        spatialPeak = get_first_available_field(All_Cells_combined(s), {'spatialPeakRate', 'SpatialPeakRate', 'spatial_peak_rate'});
        remapping = get_field_or_empty(All_Cells_combined(s), 'corr_of1_of2');
        stability = get_field_or_empty(All_Cells_combined(s), 'corr_of1_of3');
        thetaModulationScore = get_field_or_empty(All_Cells_combined(s), 'thetaMod_score');
        thetaModulationTanaka = get_field_or_empty(All_Cells_combined(s), 'mod_score_tanaka');
        intrinsicFrequency = get_field_or_empty(All_Cells_combined(s), 'intrins_frequ');
        thetaPhasePreferred = get_field_or_empty(All_Cells_combined(s), 'thetaPhase_pref');
        thetaPhaseR = get_field_or_empty(All_Cells_combined(s), 'thetaPhase_R');
        thetaPhaseP = get_field_or_empty(All_Cells_combined(s), 'thetaPhase_p');
        thetaPhaseZ = get_field_or_empty(All_Cells_combined(s), 'thetaPhase_z');
        thetaPhaseNSpikes = get_field_or_empty(All_Cells_combined(s), 'thetaPhase_nSpikes');
        thetaModulationScore1 = get_of_metric_or_empty(All_Cells_combined(s), ...
            {'thetaMod_score_of1', 'of1_thetaMod_score'}, {'thetaMod_score_sess'}, 1);
        thetaModulationScore2 = get_of_metric_or_empty(All_Cells_combined(s), ...
            {'thetaMod_score_of2', 'of2_thetaMod_score'}, {'thetaMod_score_sess'}, 2);
        thetaModulationScore3 = get_of_metric_or_empty(All_Cells_combined(s), ...
            {'thetaMod_score_of3', 'of3_thetaMod_score'}, {'thetaMod_score_sess'}, 3);
        thetaModulationTanaka1 = get_of_metric_or_empty(All_Cells_combined(s), ...
            {'mod_score_tanaka_of1', 'of1_mod_score_tanaka'}, {'mod_score_tanaka_sess'}, 1);
        thetaModulationTanaka2 = get_of_metric_or_empty(All_Cells_combined(s), ...
            {'mod_score_tanaka_of2', 'of2_mod_score_tanaka'}, {'mod_score_tanaka_sess'}, 2);
        thetaModulationTanaka3 = get_of_metric_or_empty(All_Cells_combined(s), ...
            {'mod_score_tanaka_of3', 'of3_mod_score_tanaka'}, {'mod_score_tanaka_sess'}, 3);
        intrinsicFrequency1 = get_of_metric_or_empty(All_Cells_combined(s), ...
            {'intrins_frequ_of1', 'of1_intrins_frequ'}, {'intrins_frequ_sess'}, 1);
        intrinsicFrequency2 = get_of_metric_or_empty(All_Cells_combined(s), ...
            {'intrins_frequ_of2', 'of2_intrins_frequ'}, {'intrins_frequ_sess'}, 2);
        intrinsicFrequency3 = get_of_metric_or_empty(All_Cells_combined(s), ...
            {'intrins_frequ_of3', 'of3_intrins_frequ'}, {'intrins_frequ_sess'}, 3);

        if s <= numel(saltPBySession) && ~isempty(saltPBySession{s})
            salt = saltPBySession{s};
        else
            salt = [];
        end

        nCells = max([numel(bm), numel(rel), numel(lat), numel(jit), numel(cls), numel(salt), ...
            numel(s1_ssmi), numel(s1_sfi), numel(s1_psp), numel(s1_sppr), ...
            numel(s2_ssmi), numel(s2_sfi), numel(s2_psp), numel(s2_sppr), ...
            numel(room), numel(burst), numel(speedScore), numel(peakRate), numel(avgRate), numel(avgRateOF), ...
            numel(avgRate1), numel(avgRate2), numel(avgRate3), ...
            numel(peakRate1), numel(peakRate2), numel(peakRate3), ...
            numel(temporalPeakRate1), numel(temporalPeakRate2), numel(temporalPeakRate3), ...
            numel(burst1), numel(burst2), numel(burst3), ...
            numel(pf1), numel(pf2), numel(pf3), ...
            numel(pfnum1), numel(pfnum2), numel(pfnum3), ...
            numel(zcoh1), numel(zcoh2), numel(zcoh3), ...
            numel(sparse1), numel(sparse2), numel(sparse3), ...
            numel(select1), numel(select2), numel(select3), ...
            numel(spatialInfo_A), numel(spatialInfo_B), numel(spatialInfo_A_revisit), ...
            numel(spatialInfoRate_A), numel(spatialInfoRate_B), numel(spatialInfoRate_A_revisit), ...
            numel(spatialPeak), numel(remapping), numel(stability), ...
            numel(thetaModulationScore), numel(thetaModulationTanaka), numel(intrinsicFrequency), ...
            numel(thetaPhasePreferred), numel(thetaPhaseR), numel(thetaPhaseP), ...
            numel(thetaPhaseZ), numel(thetaPhaseNSpikes), ...
            numel(thetaModulationScore1), numel(thetaModulationScore2), numel(thetaModulationScore3), ...
            numel(thetaModulationTanaka1), numel(thetaModulationTanaka2), numel(thetaModulationTanaka3), ...
            numel(intrinsicFrequency1), numel(intrinsicFrequency2), numel(intrinsicFrequency3), 0]);
        if nCells == 0
            continue
        end

        bm_pad = nan(nCells, 1);
        rel_pad = nan(nCells, 1);
        lat_pad = nan(nCells, 1);
        jit_pad = nan(nCells, 1);
        salt_pad = nan(nCells, 1);
        room_pad = nan(nCells, 1);
        burst_pad = nan(nCells, 1);
        speed_score_pad = nan(nCells, 1);
        peak_rate_pad = nan(nCells, 1);
        avg_rate_pad = nan(nCells, 1);
        avg_rate_of_pad = nan(nCells, 1);
        avg_rate1_pad = nan(nCells, 1);
        avg_rate2_pad = nan(nCells, 1);
        avg_rate3_pad = nan(nCells, 1);
        peak_rate1_pad = nan(nCells, 1);
        peak_rate2_pad = nan(nCells, 1);
        peak_rate3_pad = nan(nCells, 1);
        temporal_peak_rate1_pad = nan(nCells, 1);
        temporal_peak_rate2_pad = nan(nCells, 1);
        temporal_peak_rate3_pad = nan(nCells, 1);
        burst1_pad = nan(nCells, 1);
        burst2_pad = nan(nCells, 1);
        burst3_pad = nan(nCells, 1);
        pf1_pad = nan(nCells, 1);
        pf2_pad = nan(nCells, 1);
        pf3_pad = nan(nCells, 1);
        pfnum1_pad = nan(nCells, 1);
        pfnum2_pad = nan(nCells, 1);
        pfnum3_pad = nan(nCells, 1);
        zcoh1_pad = nan(nCells, 1);
        zcoh2_pad = nan(nCells, 1);
        zcoh3_pad = nan(nCells, 1);
        sparse1_pad = nan(nCells, 1);
        sparse2_pad = nan(nCells, 1);
        sparse3_pad = nan(nCells, 1);
        select1_pad = nan(nCells, 1);
        select2_pad = nan(nCells, 1);
        select3_pad = nan(nCells, 1);
        spatial_info_a_pad = nan(nCells, 1);
        spatial_info_b_pad = nan(nCells, 1);
        spatial_info_a_revisit_pad = nan(nCells, 1);
        spatial_info_rate_a_pad = nan(nCells, 1);
        spatial_info_rate_b_pad = nan(nCells, 1);
        spatial_info_rate_a_revisit_pad = nan(nCells, 1);
        spatial_peak_rate_pad = nan(nCells, 1);
        remapping_pad = nan(nCells, 1);
        stability_pad = nan(nCells, 1);
        theta_modulation_score_pad = nan(nCells, 1);
        theta_modulation_tanaka_pad = nan(nCells, 1);
        intrinsic_frequency_pad = nan(nCells, 1);
        theta_phase_preferred_pad = nan(nCells, 1);
        theta_phase_r_pad = nan(nCells, 1);
        theta_phase_p_pad = nan(nCells, 1);
        theta_phase_z_pad = nan(nCells, 1);
        theta_phase_n_spikes_pad = nan(nCells, 1);
        theta_modulation_score1_pad = nan(nCells, 1);
        theta_modulation_score2_pad = nan(nCells, 1);
        theta_modulation_score3_pad = nan(nCells, 1);
        theta_modulation_tanaka1_pad = nan(nCells, 1);
        theta_modulation_tanaka2_pad = nan(nCells, 1);
        theta_modulation_tanaka3_pad = nan(nCells, 1);
        intrinsic_frequency1_pad = nan(nCells, 1);
        intrinsic_frequency2_pad = nan(nCells, 1);
        intrinsic_frequency3_pad = nan(nCells, 1);
        s1_ssmi_pad = nan(nCells, 1);
        s1_sfi_pad = nan(nCells, 1);
        s1_psp_pad = nan(nCells, 1);
        s1_sppr_pad = nan(nCells, 1);
        s2_ssmi_pad = nan(nCells, 1);
        s2_sfi_pad = nan(nCells, 1);
        s2_psp_pad = nan(nCells, 1);
        s2_sppr_pad = nan(nCells, 1);
        cls_pad = repmat("unknown", nCells, 1);
        animal_pad = repmat("", nCells, 1);
        room_group_pad = repmat("unknown", nCells, 1);
        room_name_pad = repmat("unknown", nCells, 1);
        i_pad = repmat(s, nCells, 1);

        if ~isempty(bm), bm_pad(1:numel(bm)) = bm(:); end
        if ~isempty(rel), rel_pad(1:numel(rel)) = rel(:); end
        if ~isempty(lat), lat_pad(1:numel(lat)) = lat(:); end
        if ~isempty(jit), jit_pad(1:numel(jit)) = jit(:); end
        if ~isempty(salt), salt_pad(1:numel(salt)) = salt(:); end
        if ~isempty(room), room_pad(1:numel(room)) = room(:); end
        if ~isempty(burst), burst_pad(1:numel(burst)) = burst(:); end
        if ~isempty(speedScore), speed_score_pad(1:numel(speedScore)) = speedScore(:); end
        if ~isempty(peakRate), peak_rate_pad(1:numel(peakRate)) = peakRate(:); end
        if ~isempty(avgRate), avg_rate_pad(1:numel(avgRate)) = avgRate(:); end
        if ~isempty(avgRateOF), avg_rate_of_pad(1:numel(avgRateOF)) = avgRateOF(:); end
        if ~isempty(avgRate1), avg_rate1_pad(1:numel(avgRate1)) = avgRate1(:); end
        if ~isempty(avgRate2), avg_rate2_pad(1:numel(avgRate2)) = avgRate2(:); end
        if ~isempty(avgRate3), avg_rate3_pad(1:numel(avgRate3)) = avgRate3(:); end
        if ~isempty(peakRate1), peak_rate1_pad(1:numel(peakRate1)) = peakRate1(:); end
        if ~isempty(peakRate2), peak_rate2_pad(1:numel(peakRate2)) = peakRate2(:); end
        if ~isempty(peakRate3), peak_rate3_pad(1:numel(peakRate3)) = peakRate3(:); end
        if ~isempty(temporalPeakRate1), temporal_peak_rate1_pad(1:numel(temporalPeakRate1)) = temporalPeakRate1(:); end
        if ~isempty(temporalPeakRate2), temporal_peak_rate2_pad(1:numel(temporalPeakRate2)) = temporalPeakRate2(:); end
        if ~isempty(temporalPeakRate3), temporal_peak_rate3_pad(1:numel(temporalPeakRate3)) = temporalPeakRate3(:); end
        if ~isempty(burst1), burst1_pad(1:numel(burst1)) = burst1(:); end
        if ~isempty(burst2), burst2_pad(1:numel(burst2)) = burst2(:); end
        if ~isempty(burst3), burst3_pad(1:numel(burst3)) = burst3(:); end
        if ~isempty(pf1), pf1_pad(1:numel(pf1)) = pf1(:); end
        if ~isempty(pf2), pf2_pad(1:numel(pf2)) = pf2(:); end
        if ~isempty(pf3), pf3_pad(1:numel(pf3)) = pf3(:); end
        if ~isempty(pfnum1), pfnum1_pad(1:numel(pfnum1)) = pfnum1(:); end
        if ~isempty(pfnum2), pfnum2_pad(1:numel(pfnum2)) = pfnum2(:); end
        if ~isempty(pfnum3), pfnum3_pad(1:numel(pfnum3)) = pfnum3(:); end
        if ~isempty(zcoh1), zcoh1_pad(1:numel(zcoh1)) = zcoh1(:); end
        if ~isempty(zcoh2), zcoh2_pad(1:numel(zcoh2)) = zcoh2(:); end
        if ~isempty(zcoh3), zcoh3_pad(1:numel(zcoh3)) = zcoh3(:); end
        if ~isempty(sparse1), sparse1_pad(1:numel(sparse1)) = sparse1(:); end
        if ~isempty(sparse2), sparse2_pad(1:numel(sparse2)) = sparse2(:); end
        if ~isempty(sparse3), sparse3_pad(1:numel(sparse3)) = sparse3(:); end
        if ~isempty(select1), select1_pad(1:numel(select1)) = select1(:); end
        if ~isempty(select2), select2_pad(1:numel(select2)) = select2(:); end
        if ~isempty(select3), select3_pad(1:numel(select3)) = select3(:); end
        if ~isempty(spatialInfo_A), spatial_info_a_pad(1:numel(spatialInfo_A)) = spatialInfo_A(:); end
        if ~isempty(spatialInfo_B), spatial_info_b_pad(1:numel(spatialInfo_B)) = spatialInfo_B(:); end
        if ~isempty(spatialInfo_A_revisit), spatial_info_a_revisit_pad(1:numel(spatialInfo_A_revisit)) = spatialInfo_A_revisit(:); end
        if ~isempty(spatialInfoRate_A), spatial_info_rate_a_pad(1:numel(spatialInfoRate_A)) = spatialInfoRate_A(:); end
        if ~isempty(spatialInfoRate_B), spatial_info_rate_b_pad(1:numel(spatialInfoRate_B)) = spatialInfoRate_B(:); end
        if ~isempty(spatialInfoRate_A_revisit), spatial_info_rate_a_revisit_pad(1:numel(spatialInfoRate_A_revisit)) = spatialInfoRate_A_revisit(:); end
        if ~isempty(spatialPeak), spatial_peak_rate_pad(1:numel(spatialPeak)) = spatialPeak(:); end
        if ~isempty(remapping), remapping_pad(1:numel(remapping)) = remapping(:); end
        if ~isempty(stability), stability_pad(1:numel(stability)) = stability(:); end
        if ~isempty(thetaModulationScore), theta_modulation_score_pad(1:numel(thetaModulationScore)) = thetaModulationScore(:); end
        if ~isempty(thetaModulationTanaka), theta_modulation_tanaka_pad(1:numel(thetaModulationTanaka)) = thetaModulationTanaka(:); end
        if ~isempty(intrinsicFrequency), intrinsic_frequency_pad(1:numel(intrinsicFrequency)) = intrinsicFrequency(:); end
        if ~isempty(thetaPhasePreferred), theta_phase_preferred_pad(1:numel(thetaPhasePreferred)) = thetaPhasePreferred(:); end
        if ~isempty(thetaPhaseR), theta_phase_r_pad(1:numel(thetaPhaseR)) = thetaPhaseR(:); end
        if ~isempty(thetaPhaseP), theta_phase_p_pad(1:numel(thetaPhaseP)) = thetaPhaseP(:); end
        if ~isempty(thetaPhaseZ), theta_phase_z_pad(1:numel(thetaPhaseZ)) = thetaPhaseZ(:); end
        if ~isempty(thetaPhaseNSpikes), theta_phase_n_spikes_pad(1:numel(thetaPhaseNSpikes)) = thetaPhaseNSpikes(:); end
        if ~isempty(thetaModulationScore1), theta_modulation_score1_pad(1:numel(thetaModulationScore1)) = thetaModulationScore1(:); end
        if ~isempty(thetaModulationScore2), theta_modulation_score2_pad(1:numel(thetaModulationScore2)) = thetaModulationScore2(:); end
        if ~isempty(thetaModulationScore3), theta_modulation_score3_pad(1:numel(thetaModulationScore3)) = thetaModulationScore3(:); end
        if ~isempty(thetaModulationTanaka1), theta_modulation_tanaka1_pad(1:numel(thetaModulationTanaka1)) = thetaModulationTanaka1(:); end
        if ~isempty(thetaModulationTanaka2), theta_modulation_tanaka2_pad(1:numel(thetaModulationTanaka2)) = thetaModulationTanaka2(:); end
        if ~isempty(thetaModulationTanaka3), theta_modulation_tanaka3_pad(1:numel(thetaModulationTanaka3)) = thetaModulationTanaka3(:); end
        if ~isempty(intrinsicFrequency1), intrinsic_frequency1_pad(1:numel(intrinsicFrequency1)) = intrinsicFrequency1(:); end
        if ~isempty(intrinsicFrequency2), intrinsic_frequency2_pad(1:numel(intrinsicFrequency2)) = intrinsicFrequency2(:); end
        if ~isempty(intrinsicFrequency3), intrinsic_frequency3_pad(1:numel(intrinsicFrequency3)) = intrinsicFrequency3(:); end
        if ~isempty(s1_ssmi), s1_ssmi_pad(1:numel(s1_ssmi)) = s1_ssmi(:); end
        if ~isempty(s1_sfi), s1_sfi_pad(1:numel(s1_sfi)) = s1_sfi(:); end
        if ~isempty(s1_psp), s1_psp_pad(1:numel(s1_psp)) = s1_psp(:); end
        if ~isempty(s1_sppr), s1_sppr_pad(1:numel(s1_sppr)) = s1_sppr(:); end
        if ~isempty(s2_ssmi), s2_ssmi_pad(1:numel(s2_ssmi)) = s2_ssmi(:); end
        if ~isempty(s2_sfi), s2_sfi_pad(1:numel(s2_sfi)) = s2_sfi(:); end
        if ~isempty(s2_psp), s2_psp_pad(1:numel(s2_psp)) = s2_psp(:); end
        if ~isempty(s2_sppr), s2_sppr_pad(1:numel(s2_sppr)) = s2_sppr(:); end
        if ~isempty(iNumber)
            iValues = double(iNumber(:));
            i_pad(1:min(numel(iValues), nCells)) = iValues(1:min(numel(iValues), nCells));
        end
        if ~isempty(animal)
            animal_pad(:) = normalize_scalar_identifier(animal);
        end

        if ~isempty(cls)
            nCls = min(numel(cls), nCells);
            for k = 1:nCls
                cls_pad(k) = map_cell_class(cls(k), principal_codes, interneuron_codes);
            end
        end

        room_group_pad(room_pad == 15) = "main";
        room_group_pad(room_pad == 19 | room_pad == 20) = "control";
        room_name_pad = room_id_display_name(room_pad);
        pf_mean_pad = row_nanmean([pf1_pad pf2_pad pf3_pad]);
        pfnum_mean_pad = row_nanmean([pfnum1_pad pfnum2_pad pfnum3_pad]);
        avg_rate1_export_pad = row_nanmean_with_fallback(avg_rate1_pad, avg_rate_of_pad);
        avg_rate2_export_pad = row_nanmean_with_fallback(avg_rate2_pad, avg_rate_of_pad);
        avg_rate3_export_pad = row_nanmean_with_fallback(avg_rate3_pad, avg_rate_of_pad);
        avg_rate_of_mean_pad = row_nanmean_with_fallback([avg_rate1_pad avg_rate2_pad avg_rate3_pad], avg_rate_of_pad);
        temporal_peak_rate1_export_pad = row_nanmean_with_fallback(temporal_peak_rate1_pad, peak_rate_pad);
        temporal_peak_rate2_export_pad = row_nanmean_with_fallback(temporal_peak_rate2_pad, peak_rate_pad);
        temporal_peak_rate3_export_pad = row_nanmean_with_fallback(temporal_peak_rate3_pad, peak_rate_pad);
        spatial_peak_rate_of_max_pad = row_nanmax_with_fallback([peak_rate1_pad peak_rate2_pad peak_rate3_pad], spatial_peak_rate_pad);
        peak_rate_of_max_pad = row_nanmax_with_fallback([peak_rate1_pad peak_rate2_pad peak_rate3_pad], peak_rate_pad);
        burst_index_of_mean_pad = row_nanmean_with_fallback([burst1_pad burst2_pad burst3_pad], burst_pad);
        theta_modulation_score_of_mean_pad = row_nanmean_with_fallback( ...
            [theta_modulation_score1_pad theta_modulation_score2_pad theta_modulation_score3_pad], ...
            theta_modulation_score_pad);
        theta_modulation_tanaka_of_mean_pad = row_nanmean_with_fallback( ...
            [theta_modulation_tanaka1_pad theta_modulation_tanaka2_pad theta_modulation_tanaka3_pad], ...
            theta_modulation_tanaka_pad);
        intrinsic_frequency_of_mean_pad = row_nanmean_with_fallback( ...
            [intrinsic_frequency1_pad intrinsic_frequency2_pad intrinsic_frequency3_pad], ...
            intrinsic_frequency_pad);

        SessionIndex = [SessionIndex; repmat(s, nCells, 1)]; %#ok<AGROW>
        INumber = [INumber; i_pad]; %#ok<AGROW>
        Cell = [Cell; (1:nCells)']; %#ok<AGROW>
        CellClass = [CellClass; cls_pad]; %#ok<AGROW>
        AnimalID = [AnimalID; animal_pad]; %#ok<AGROW>
        RoomID = [RoomID; room_pad]; %#ok<AGROW>
        RoomName = [RoomName; room_name_pad]; %#ok<AGROW>
        RoomGroup = [RoomGroup; room_group_pad]; %#ok<AGROW>
        BaselineSD = [BaselineSD; bm_pad]; %#ok<AGROW>
        Reliability = [Reliability; rel_pad]; %#ok<AGROW>
        Latency_ms = [Latency_ms; lat_pad]; %#ok<AGROW>
        Jitter_ms = [Jitter_ms; jit_pad]; %#ok<AGROW>
        SALTPValue = [SALTPValue; salt_pad]; %#ok<AGROW>
        S1_SSMI = [S1_SSMI; s1_ssmi_pad]; %#ok<AGROW>
        S1_SFI = [S1_SFI; s1_sfi_pad]; %#ok<AGROW>
        S1_PSP = [S1_PSP; s1_psp_pad]; %#ok<AGROW>
        S1_SpPR = [S1_SpPR; s1_sppr_pad]; %#ok<AGROW>
        S2_SSMI = [S2_SSMI; s2_ssmi_pad]; %#ok<AGROW>
        S2_SFI = [S2_SFI; s2_sfi_pad]; %#ok<AGROW>
        S2_PSP = [S2_PSP; s2_psp_pad]; %#ok<AGROW>
        S2_SpPR = [S2_SpPR; s2_sppr_pad]; %#ok<AGROW>
        BurstIndex = [BurstIndex; burst_pad]; %#ok<AGROW>
        SpeedScore = [SpeedScore; speed_score_pad]; %#ok<AGROW>
        PeakRate = [PeakRate; peak_rate_pad]; %#ok<AGROW>
        TemporalPeakRate = [TemporalPeakRate; peak_rate_pad]; %#ok<AGROW>
        AverageRate = [AverageRate; avg_rate_pad]; %#ok<AGROW>
        AverageRate_OF1 = [AverageRate_OF1; avg_rate1_export_pad]; %#ok<AGROW>
        AverageRate_OF2 = [AverageRate_OF2; avg_rate2_export_pad]; %#ok<AGROW>
        AverageRate_OF3 = [AverageRate_OF3; avg_rate3_export_pad]; %#ok<AGROW>
        PlaceFieldSize = [PlaceFieldSize; pf_mean_pad]; %#ok<AGROW>
        SpatialInfo_A = [SpatialInfo_A; spatial_info_a_pad]; %#ok<AGROW>
        SpatialPeakRate = [SpatialPeakRate; spatial_peak_rate_of_max_pad]; %#ok<AGROW>
        SpatialPeakRate_OF1 = [SpatialPeakRate_OF1; peak_rate1_pad]; %#ok<AGROW>
        SpatialPeakRate_OF2 = [SpatialPeakRate_OF2; peak_rate2_pad]; %#ok<AGROW>
        SpatialPeakRate_OF3 = [SpatialPeakRate_OF3; peak_rate3_pad]; %#ok<AGROW>
        PlaceFieldSize_OF1 = [PlaceFieldSize_OF1; pf1_pad]; %#ok<AGROW>
        PlaceFieldSize_OF2 = [PlaceFieldSize_OF2; pf2_pad]; %#ok<AGROW>
        PlaceFieldSize_OF3 = [PlaceFieldSize_OF3; pf3_pad]; %#ok<AGROW>
        PlaceFieldNumber = [PlaceFieldNumber; pfnum_mean_pad]; %#ok<AGROW>
        PlaceFieldNumber_OF1 = [PlaceFieldNumber_OF1; pfnum1_pad]; %#ok<AGROW>
        PlaceFieldNumber_OF2 = [PlaceFieldNumber_OF2; pfnum2_pad]; %#ok<AGROW>
        PlaceFieldNumber_OF3 = [PlaceFieldNumber_OF3; pfnum3_pad]; %#ok<AGROW>
        ZCoherence_OF1 = [ZCoherence_OF1; zcoh1_pad]; %#ok<AGROW>
        ZCoherence_OF2 = [ZCoherence_OF2; zcoh2_pad]; %#ok<AGROW>
        ZCoherence_OF3 = [ZCoherence_OF3; zcoh3_pad]; %#ok<AGROW>
        Sparseness_OF1 = [Sparseness_OF1; sparse1_pad]; %#ok<AGROW>
        Sparseness_OF2 = [Sparseness_OF2; sparse2_pad]; %#ok<AGROW>
        Sparseness_OF3 = [Sparseness_OF3; sparse3_pad]; %#ok<AGROW>
        Selectivity_OF1 = [Selectivity_OF1; select1_pad]; %#ok<AGROW>
        Selectivity_OF2 = [Selectivity_OF2; select2_pad]; %#ok<AGROW>
        Selectivity_OF3 = [Selectivity_OF3; select3_pad]; %#ok<AGROW>
        SpatialInfo_OF1 = [SpatialInfo_OF1; spatial_info_a_pad]; %#ok<AGROW>
        SpatialInfo_OF2 = [SpatialInfo_OF2; spatial_info_b_pad]; %#ok<AGROW>
        SpatialInfo_OF3 = [SpatialInfo_OF3; spatial_info_a_revisit_pad]; %#ok<AGROW>
        SpatialInfoRate_OF1 = [SpatialInfoRate_OF1; spatial_info_rate_a_pad]; %#ok<AGROW>
        SpatialInfoRate_OF2 = [SpatialInfoRate_OF2; spatial_info_rate_b_pad]; %#ok<AGROW>
        SpatialInfoRate_OF3 = [SpatialInfoRate_OF3; spatial_info_rate_a_revisit_pad]; %#ok<AGROW>
        Remapping = [Remapping; remapping_pad]; %#ok<AGROW>
        Stability = [Stability; stability_pad]; %#ok<AGROW>
        ThetaModulationScore = [ThetaModulationScore; theta_modulation_score_pad]; %#ok<AGROW>
        ThetaModulationTanaka = [ThetaModulationTanaka; theta_modulation_tanaka_pad]; %#ok<AGROW>
        IntrinsicFrequency = [IntrinsicFrequency; intrinsic_frequency_pad]; %#ok<AGROW>
        ThetaPhasePreferred = [ThetaPhasePreferred; theta_phase_preferred_pad]; %#ok<AGROW>
        ThetaPhaseR = [ThetaPhaseR; theta_phase_r_pad]; %#ok<AGROW>
        ThetaPhaseP = [ThetaPhaseP; theta_phase_p_pad]; %#ok<AGROW>
        ThetaPhaseZ = [ThetaPhaseZ; theta_phase_z_pad]; %#ok<AGROW>
        ThetaPhaseNSpikes = [ThetaPhaseNSpikes; theta_phase_n_spikes_pad]; %#ok<AGROW>
        AverageRate_OFMean = [AverageRate_OFMean; avg_rate_of_mean_pad]; %#ok<AGROW>
        PeakRate_OFMax = [PeakRate_OFMax; peak_rate_of_max_pad]; %#ok<AGROW>
        SpatialPeakRate_OFMax = [SpatialPeakRate_OFMax; spatial_peak_rate_of_max_pad]; %#ok<AGROW>
        TemporalPeakRate_OF1 = [TemporalPeakRate_OF1; temporal_peak_rate1_export_pad]; %#ok<AGROW>
        TemporalPeakRate_OF2 = [TemporalPeakRate_OF2; temporal_peak_rate2_export_pad]; %#ok<AGROW>
        TemporalPeakRate_OF3 = [TemporalPeakRate_OF3; temporal_peak_rate3_export_pad]; %#ok<AGROW>
        BurstIndex_OFMean = [BurstIndex_OFMean; burst_index_of_mean_pad]; %#ok<AGROW>
        ThetaModulationScore_OFMean = [ThetaModulationScore_OFMean; theta_modulation_score_of_mean_pad]; %#ok<AGROW>
        ThetaModulationTanaka_OFMean = [ThetaModulationTanaka_OFMean; theta_modulation_tanaka_of_mean_pad]; %#ok<AGROW>
        IntrinsicFrequency_OFMean = [IntrinsicFrequency_OFMean; intrinsic_frequency_of_mean_pad]; %#ok<AGROW>
    end

    T = table(SessionIndex, INumber, Cell, CellClass, AnimalID, RoomID, RoomName, RoomGroup, BaselineSD, Reliability, ...
        Latency_ms, Jitter_ms, SALTPValue, S1_SSMI, S1_SFI, S1_PSP, S1_SpPR, ...
        S2_SSMI, S2_SFI, S2_PSP, S2_SpPR, ...
        BurstIndex, SpeedScore, PeakRate, TemporalPeakRate, ...
        AverageRate, AverageRate_OF1, AverageRate_OF2, AverageRate_OF3, ...
        PlaceFieldSize, SpatialInfo_A, SpatialPeakRate, ...
        SpatialPeakRate_OF1, SpatialPeakRate_OF2, SpatialPeakRate_OF3, ...
        PlaceFieldSize_OF1, PlaceFieldSize_OF2, PlaceFieldSize_OF3, ...
        PlaceFieldNumber, ...
        PlaceFieldNumber_OF1, PlaceFieldNumber_OF2, PlaceFieldNumber_OF3, ...
        ZCoherence_OF1, ZCoherence_OF2, ZCoherence_OF3, ...
        Sparseness_OF1, Sparseness_OF2, Sparseness_OF3, ...
        Selectivity_OF1, Selectivity_OF2, Selectivity_OF3, ...
        SpatialInfo_OF1, SpatialInfo_OF2, SpatialInfo_OF3, ...
        SpatialInfoRate_OF1, SpatialInfoRate_OF2, SpatialInfoRate_OF3, ...
        Remapping, Stability, ThetaModulationScore, ThetaModulationTanaka, IntrinsicFrequency, ...
        ThetaPhasePreferred, ThetaPhaseR, ThetaPhaseP, ThetaPhaseZ, ThetaPhaseNSpikes, ...
        AverageRate_OFMean, PeakRate_OFMax, SpatialPeakRate_OFMax, ...
        TemporalPeakRate_OF1, TemporalPeakRate_OF2, TemporalPeakRate_OF3, ...
        BurstIndex_OFMean, ...
        ThetaModulationScore_OFMean, ThetaModulationTanaka_OFMean, IntrinsicFrequency_OFMean);
end

function [Tin, Thresholds] = run_salt_two_step_classification(Tin, default_p_cutoff, min_reliability_for_tagged, stageB_settings)

    n = height(Tin);
    knownClassMask = Tin.CellClass == "principal" | Tin.CellClass == "interneuron";
    validP = ~isnan(Tin.SALTPValue) & knownClassMask;
    reliableEnoughMask = ~isnan(Tin.Reliability) & Tin.Reliability >= min_reliability_for_tagged;
    taggedMask = validP & reliableEnoughMask & Tin.SALTPValue < default_p_cutoff;
    untaggedMask = validP & ~taggedMask;
    lowReliabilityTaggedCandidates = validP & Tin.SALTPValue < default_p_cutoff & ~reliableEnoughMask;

    Tin.StageA_Label = repmat("unclassified", n, 1);
    Tin.StageA_Label(untaggedMask) = "untagged";
    Tin.StageA_Label(taggedMask) = "tagged";

    Tin.StageB_Label = repmat("indirect", n, 1);
    Tin.FinalLabel = repmat("unclassified", n, 1);
    Tin.P_direct_given_tagged = nan(n, 1);
    Tin.DirectScore = nan(n, 1);

    Thresholds = struct();
    Thresholds.StageA = struct();
    Thresholds.StageA.PValueCutoff = default_p_cutoff;
    Thresholds.StageA.MinReliabilityForTagged = min_reliability_for_tagged;
    Thresholds.StageA.Diagnostics = struct();
    Thresholds.StageA.Diagnostics.NumValidPValue = nnz(validP);
    Thresholds.StageA.Diagnostics.NumReliableEnough = nnz(validP & reliableEnoughMask);
    Thresholds.StageA.Diagnostics.NumBelowReliabilityCutoff = nnz(validP & ~reliableEnoughMask);
    Thresholds.StageA.Diagnostics.NumPValuePositiveButReliabilityRejected = nnz(lowReliabilityTaggedCandidates);
    Thresholds.StageA.Diagnostics.NumTagged = nnz(taggedMask);
    Thresholds.StageA.Diagnostics.NumUntagged = nnz(untaggedMask);
    Thresholds.StageA.Diagnostics.FracBelow0_05 = mean(Tin.SALTPValue(validP) < 0.05, 'omitnan');
    Thresholds.StageA.Diagnostics.FracBelow0_01 = mean(Tin.SALTPValue(validP) < 0.01, 'omitnan');
    Thresholds.StageA.Diagnostics.FracBelowCutoff = mean(Tin.SALTPValue(validP) < default_p_cutoff, 'omitnan');
    Thresholds.StageA.Diagnostics.FracTaggedAfterReliabilityGate = mean(taggedMask(validP), 'omitnan');

    [Tin, Thresholds.StageB] = run_salt_latency_stage_b(Tin, 'all cells', stageB_settings);

    Tin.FinalLabel(untaggedMask) = "untagged";
    Tin.FinalLabel(taggedMask) = "indirect";
    Tin.FinalLabel(Tin.StageB_Label == "direct") = "direct";

    fprintf('\n[SALT] Step A diagnostics\n');
    fprintf('  SALT p-value cutoff: %.4g\n', default_p_cutoff);
    fprintf('  Reliability cutoff for tagging: %.3f\n', min_reliability_for_tagged);
    fprintf('  Valid p-values: %d\n', nnz(validP));
    fprintf('  Rejected by reliability gate: %d\n', nnz(lowReliabilityTaggedCandidates));
    fprintf('  Tagged: %d\n', nnz(taggedMask));
    fprintf('  Untagged: %d\n', nnz(untaggedMask));
end

function add_analysis_paths(custom_settings)

    codeFolder = fileparts(mfilename('fullpath'));
    matlabCodeFolder = fullfile(codeFolder, '..');
    if exist(matlabCodeFolder, 'dir')
        addpath(genpath(matlabCodeFolder));
    end

    extraPaths = get_override_value(custom_settings, 'additionalPaths', {});
    if ischar(extraPaths) || isstring(extraPaths)
        extraPaths = cellstr(string(extraPaths));
    end
    for iPath = 1:numel(extraPaths)
        thisPath = char(extraPaths{iPath});
        if exist(thisPath, 'dir')
            addpath(genpath(thisPath));
        end
    end
end

function [cacheOk, cachedMetrics, message] = try_load_cached_salt_metrics( ...
    saltMetricPath, nCells, response_window_ms, delay_ms, guard_ms, ...
    latency_bin_ms, baseline_bootstrap_sets, null_pair_draws, ...
    min_baseline_windows, analyze_dir)

    cacheOk = false;
    message = '';
    cachedMetrics = struct( ...
        'salt_p_value', nan(nCells, 1), ...
        'salt_z_value', nan(nCells, 1), ...
        'salt_statistic', nan(nCells, 1));

    try
        S = load(saltMetricPath, 'sessionMetrics');
    catch ME
        message = sprintf('Could not load cached SALT metrics: %s', ME.message);
        return
    end

    if ~isfield(S, 'sessionMetrics')
        message = 'Cached SALT file missing sessionMetrics. Recomputing.';
        return
    end

    sessionMetrics = S.sessionMetrics;
    if ~isfield(sessionMetrics, 'salt_p_value')
        message = 'Cached SALT file missing p-values. Recomputing.';
        return
    end

    pVals = sessionMetrics.salt_p_value(:);
    if numel(pVals) ~= nCells
        message = sprintf('Cached SALT metrics have %d cells but tList has %d. Recomputing.', numel(pVals), nCells);
        return
    end

    paramsOk = true;
    if isfield(sessionMetrics, 'parameters')
        paramsOk = salt_cache_parameters_match(sessionMetrics.parameters, ...
            response_window_ms, delay_ms, guard_ms, latency_bin_ms, ...
            baseline_bootstrap_sets, null_pair_draws, min_baseline_windows, analyze_dir);
    end

    if ~paramsOk
        message = 'Cached SALT metrics were generated with different parameters. Recomputing.';
        return
    end

    cachedMetrics.salt_p_value = pVals;
    if isfield(sessionMetrics, 'salt_z_value') && numel(sessionMetrics.salt_z_value) == nCells
        cachedMetrics.salt_z_value = sessionMetrics.salt_z_value(:);
    end
    if isfield(sessionMetrics, 'salt_statistic') && numel(sessionMetrics.salt_statistic) == nCells
        cachedMetrics.salt_statistic = sessionMetrics.salt_statistic(:);
    end

    cacheOk = true;
    message = sprintf('Loaded cached SALT metrics from %s', saltMetricPath);
end

function tt_files = read_tt_list(tListPath)

    tt_files = {};
    fid = fopen(tListPath);
    if fid == -1
        return
    end

    cleaner = onCleanup(@() fclose(fid)); %#ok<NASGU>

    cell_no = 0;
    while true
        tline = fgetl(fid);
        if ~ischar(tline)
            break
        end
        cell_no = cell_no + 1;
        tt_files{cell_no} = tline; %#ok<AGROW>
    end
end

function txt = normalize_text_field(value)

    txt = '';

    if isempty(value)
        return
    end

    if isstring(value)
        value = value(1);
        if strlength(value) == 0
            return
        end
        txt = char(value);
        return
    end

    if ischar(value)
        txt = value;
        return
    end

    if iscell(value)
        if isempty(value)
            return
        end
        txt = normalize_text_field(value{1});
    end
end

function windows = build_baseline_windows(ttlOn, ttlOff, response_window_s, guard_s)

    windows = zeros(0, 2);
    if numel(ttlOn) < 2
        return
    end

    for p = 1:(numel(ttlOn) - 1)
        intervalStart = ttlOff(p) + guard_s;
        intervalEnd = ttlOn(p + 1) - guard_s;
        intervalLength = intervalEnd - intervalStart;

        if intervalLength < response_window_s
            continue
        end

        nWindows = floor(intervalLength / response_window_s);
        for w = 1:nWindows
            wStart = intervalStart + (w - 1) * response_window_s;
            wEnd = wStart + response_window_s;
            windows(end + 1, :) = [wStart, wEnd]; %#ok<AGROW>
        end
    end
end

function outDir = resolve_output_dir(preferredOutDir, stageB_settings)

    if nargin < 1 || isempty(preferredOutDir)
        preferredOutDir = fullfile(pwd, 'OptoMetricComparison_SALT');
    end

    outDir = preferredOutDir;
    parentDir = fileparts(preferredOutDir);

    if ~(isempty(parentDir) || exist(parentDir, 'dir'))
        outDir = fullfile(pwd, 'OptoMetricComparison_SALT');
    end

    outDir = fullfile(outDir, stage_b_export_folder_name(stageB_settings));
end

function allCellsPath = resolve_all_cells_path()

    codeFolder = fileparts(mfilename('fullpath'));
    repoRoot = fullfile(codeFolder, '..', '..');
    candidatePaths = { ...
        fullfile(repoRoot, 'All_Cells_combined.mat'), ...
        fullfile(repoRoot, 'Data', 'All_Cells_combined.mat'), ...
        fullfile(pwd, 'All_Cells_combined.mat')};

    allCellsPath = '';
    for i = 1:numel(candidatePaths)
        if exist(candidatePaths{i}, 'file')
            allCellsPath = candidatePaths{i};
            return
        end
    end

    error('Could not find All_Cells_combined.mat. Pass custom_settings.allCellsPath or place it in the repository root/Data folder.')
end

function sessionInfoPath = resolve_session_info_path()

    codeFolder = fileparts(mfilename('fullpath'));
    repoRoot = fullfile(codeFolder, '..', '..');
    candidatePaths = { ...
        fullfile(repoRoot, 'sessionInfo.mat'), ...
        fullfile(repoRoot, 'Data', 'sessionInfo.mat'), ...
        fullfile(repoRoot, 'Analysis_scripts', 'DataOrganization', 'sessionInfo.mat'), ...
        fullfile(pwd, 'sessionInfo.mat')};

    sessionInfoPath = '';
    for i = 1:numel(candidatePaths)
        if exist(candidatePaths{i}, 'file')
            sessionInfoPath = candidatePaths{i};
            return
        end
    end

    error('Could not find sessionInfo.mat. Pass custom_settings.sessionInfoPath or place it in the repository root/Data folder.')
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

function stageB_settings = ensure_stage_b_settings_defaults(stageB_settings)

    if nargin < 1 || isempty(stageB_settings)
        stageB_settings = struct();
    end

    if ~isfield(stageB_settings, 'fixed_latency_threshold_ms') || isempty(stageB_settings.fixed_latency_threshold_ms)
        stageB_settings.fixed_latency_threshold_ms = 7;
    end

    if ~isfield(stageB_settings, 'principal_latency_hist_bin_method') || isempty(stageB_settings.principal_latency_hist_bin_method)
        stageB_settings.principal_latency_hist_bin_method = 'integers';
    end
end

function description = describe_stage_b_export_variant(stageB_settings)

    stageB_settings = ensure_stage_b_settings_defaults(stageB_settings);
    description = sprintf('fixed latency threshold at %s ms; latency >= threshold is indirect', ...
        sanitize_threshold_text(stageB_settings.fixed_latency_threshold_ms));
end

function run_salt_procedure_figures(T, SaltTwoStep, ExportSettings, outDir)

    if isempty(T) || height(T) == 0
        return
    end

    default_salt_cutoff = ExportSettings.default_salt_cutoff;
    stageB_settings = ensure_stage_b_settings_defaults(ExportSettings.stageB_settings);
    figureRoot = fullfile(outDir, 'SALT_procedure_figures');
    if ~exist(figureRoot, 'dir')
        mkdir(figureRoot);
    end

    classDefs = { ...
        'principal', T.CellClass == "principal"; ...
        'interneuron', T.CellClass == "interneuron"; ...
        'all_classes', true(height(T), 1)};

    for iClass = 1:size(classDefs, 1)
        className = classDefs{iClass, 1};
        classMask = classDefs{iClass, 2};
        if nnz(classMask) == 0
            continue
        end

        classOutDir = fullfile(figureRoot, sanitize_filename(className));
        if ~exist(classOutDir, 'dir')
            mkdir(classOutDir);
        end

        make_pvalue_distribution_plot(T(classMask,:), className, classOutDir, default_salt_cutoff);
        make_tagged_latency_histogram(T(classMask,:), className, classOutDir, SaltTwoStep.StageB, stageB_settings);
    end

    export_salt_manuscript_figure_outputs(T, SaltTwoStep.StageB, outDir, ExportSettings);
end

function last_updated = mark_field_updated(last_updated, session_index, field_name, update_timestamp, update_script_name)

    last_updated(session_index).(field_name) = struct( ...
        'update_timestamp', update_timestamp, ...
        'update_script', update_script_name);
end

function h = category_histogram(cats, nCategories)

    counts = accumarray(cats(:), 1, [nCategories 1], @sum, 0);
    if sum(counts) == 0
        h = ones(1, nCategories) / nCategories;
    else
        h = (counts ./ sum(counts))';
    end
end

function d = jensen_shannon_divergence(p, q)

    p = p(:)';
    q = q(:)';
    p = p / sum(p);
    q = q / sum(q);
    m = 0.5 * (p + q);

    d = 0.5 * kl_divergence(p, m) + 0.5 * kl_divergence(q, m);
end

function y = row_nanmean(X)

    validCounts = sum(~isnan(X), 2);
    X0 = X;
    X0(isnan(X0)) = 0;

    y = nan(size(X, 1), 1);
    validRows = validCounts > 0;
    y(validRows) = sum(X0(validRows, :), 2) ./ validCounts(validRows);
end

function y = row_nanmean_with_fallback(X, fallback)

    y = row_nanmean(X);
    if nargin < 2 || isempty(fallback)
        return
    end

    fallback = fallback(:);
    missingRows = ~any(isfinite(X), 2);
    y(missingRows) = fallback(missingRows);
end

function y = row_nanmax_with_fallback(X, fallback)

    y = nan(size(X, 1), 1);
    validRows = any(isfinite(X), 2);
    X0 = X;
    X0(~isfinite(X0)) = -Inf;
    y(validRows) = max(X0(validRows, :), [], 2);

    if nargin < 2 || isempty(fallback)
        return
    end

    fallback = fallback(:);
    y(~validRows) = fallback(~validRows);
end

function x = get_field_or_empty(S, fieldname)
    if isfield(S, fieldname)
        x = S.(fieldname);
        if isempty(x)
            x = [];
        end
    else
        x = [];
    end
end

function x = get_first_available_field(S, fieldnames)

    x = [];
    for iField = 1:numel(fieldnames)
        x = get_field_or_empty(S, fieldnames{iField});
        if ~isempty(x)
            return
        end
    end
end

function x = get_temporal_peak_rate_of(S, ofIdx)

    scalarFieldNames = { ...
        sprintf('peak_fr_rate_of%d', ofIdx), ...
        sprintf('of%d_peak_fr_rate', ofIdx), ...
        sprintf('temporalPeakRate_of%d', ofIdx), ...
        sprintf('TemporalPeakRate_OF%d', ofIdx), ...
        sprintf('temporal_peak_rate_of%d', ofIdx), ...
        sprintf('of%d_temporalPeakRate', ofIdx), ...
        sprintf('of%d_peak_Rate', ofIdx)};
    matrixFieldNames = { ...
        'fields_peak_fr_rate', ...
        'of_peak_rate', ...
        'of_peak_fr_rate', ...
        'open_fields_peak_fr_rate', ...
        'peak_fr_rate_sess', ...
        'temporalPeakRate_sess', ...
        'temporal_peak_rate_sess', ...
        'peak_rate_by_of'};

    x = get_of_metric_or_empty(S, scalarFieldNames, matrixFieldNames, ofIdx);
end

function x = get_of_metric_or_empty(S, scalarFieldNames, matrixFieldNames, ofIdx)

    x = get_first_available_field(S, scalarFieldNames);
    if ~isempty(x)
        return
    end

    x = [];
    for iField = 1:numel(matrixFieldNames)
        values = get_field_or_empty(S, matrixFieldNames{iField});
        if isempty(values)
            continue
        end

        if iscell(values)
            if numel(values) >= ofIdx && isnumeric(values{ofIdx}) && ~isempty(values{ofIdx})
                x = values{ofIdx}(:);
                return
            end
            continue
        end

        if ndims(values) > 2 || isvector(values)
            continue
        end

        if size(values, 2) >= ofIdx
            x = values(:, ofIdx);
            return
        end

        if size(values, 1) >= ofIdx
            x = values(ofIdx, :).';
            return
        end
    end
end

function id = normalize_scalar_identifier(value)

    id = "";

    if isempty(value)
        return
    end

    if iscell(value)
        value = value{1};
        if isempty(value)
            return
        end
    end

    if isstring(value)
        value = value(1);
        id = strtrim(value);
        return
    end

    if ischar(value)
        id = string(strtrim(value));
        return
    end

    if isnumeric(value) || islogical(value)
        value = value(1);
        if isfinite(value)
            id = string(value);
        end
        return
    end

    try
        id = string(value);
    catch
        id = "";
    end
end

function cls = map_cell_class(val, principal_codes, interneuron_codes)

    cls = "unknown";

    if isempty(val)
        return
    end

    if iscell(val)
        if isempty(val{1})
            return
        end
        val = val{1};
    end

    if isnumeric(val) || islogical(val)
        if any(val == principal_codes)
            cls = "principal";
        elseif any(val == interneuron_codes)
            cls = "interneuron";
        end
        return
    end

    val = string(val);
    v = lower(strtrim(val));

    if contains(v, "principal") || contains(v, "pyr") || strcmp(v, "pc")
        cls = "principal";
    elseif contains(v, "interneuron") || strcmp(v, "in")
        cls = "interneuron";
    end
end

function [Tin, StageB] = run_salt_latency_stage_b(Tin, className, stageB_settings)

    stageB_settings = ensure_stage_b_settings_defaults(stageB_settings);
    taggedMask = Tin.StageA_Label == "tagged";
    validB = taggedMask & ~isnan(Tin.Latency_ms);
    threshold_ms = stageB_settings.fixed_latency_threshold_ms;

    Tin.StageB_Label(taggedMask) = "direct";
    Tin.P_direct_given_tagged(taggedMask) = NaN;

    StageB = initialize_salt_stage_b_struct('fixed_latency_threshold');
    StageB.LatencyThreshold_ms = threshold_ms;
    StageB.UsedLateLatencySplit = true;
    StageB.Diagnostics.NumTagged = nnz(taggedMask);
    StageB.Diagnostics.NumWithLatency = nnz(validB);
    StageB.Diagnostics.ClassificationMode = "fixed_latency_threshold";
    StageB.Diagnostics.SharedAcrossClasses = true;
    StageB.Diagnostics.FixedLatencyThreshold_ms = threshold_ms;
    StageB.Diagnostics.FixedThresholdOverrideApplied = true;
    StageB.Diagnostics.IndirectRule = sprintf('Latency_ms >= %.3g', threshold_ms);

    if nnz(validB) > 0
        validIdxB = find(validB);
        latValid = Tin.Latency_ms(validB);
        indirectByRule = latValid >= threshold_ms;
        directByRule = latValid < threshold_ms;

        Tin.StageB_Label(validIdxB(indirectByRule)) = "indirect";
        Tin.P_direct_given_tagged(validIdxB(directByRule)) = 1;
        Tin.P_direct_given_tagged(validIdxB(indirectByRule)) = 0;
        Tin.DirectScore(validB) = threshold_ms - latValid;
    else
        warning('%s: no tagged cells with latency. Keeping tagged cells direct.', className)
    end

    StageB = populate_salt_stage_b_raw_medians(StageB, Tin, taggedMask);
    StageB.FixedThresholdSeparation = summarize_fixed_threshold_separation(Tin, threshold_ms);

    fprintf('\n[%s] SALT Step B diagnostics\n', className);
    fprintf('  Method: fixed latency threshold\n');
    fprintf('  Tagged: %d\n', StageB.Diagnostics.NumTagged);
    fprintf('  Tagged with latency: %d\n', StageB.Diagnostics.NumWithLatency);
    fprintf('  Direct latency rule: Latency < %.4f ms\n', threshold_ms);
    fprintf('  Indirect latency rule: Latency >= %.4f ms\n', threshold_ms);
    fprintf('  Direct: %d\n', nnz(taggedMask & (Tin.StageB_Label == "direct")));
    fprintf('  Indirect: %d\n', nnz(taggedMask & (Tin.StageB_Label == "indirect")));
    fprintf('  Fixed-threshold Latency Ashman D: %.3f\n', StageB.FixedThresholdSeparation.LatencyAshmanD);
    fprintf('  Fixed-threshold mean latency gap: %.3f ms\n', StageB.FixedThresholdSeparation.LatencyMeanGap_ms);
end

function tf = salt_cache_parameters_match(parameters, response_window_ms, delay_ms, guard_ms, ...
    latency_bin_ms, baseline_bootstrap_sets, null_pair_draws, min_baseline_windows, analyze_dir)

    tf = true;

    expected = { ...
        'response_window_ms', response_window_ms; ...
        'delay_ms', delay_ms; ...
        'guard_ms', guard_ms; ...
        'latency_bin_ms', latency_bin_ms; ...
        'baseline_bootstrap_sets', baseline_bootstrap_sets; ...
        'null_pair_draws', null_pair_draws; ...
        'min_baseline_windows', min_baseline_windows};

    for i = 1:size(expected, 1)
        fieldName = expected{i,1};
        fieldValue = expected{i,2};
        if ~isfield(parameters, fieldName) || ~isequal(parameters.(fieldName), fieldValue)
            tf = false;
            return
        end
    end

    if ~isfield(parameters, 'analyze_dir')
        return
    end

    tf = strcmp(normalize_text_field(parameters.analyze_dir), analyze_dir);
end

function folderName = stage_b_export_folder_name(stageB_settings)

    stageB_settings = ensure_stage_b_settings_defaults(stageB_settings);
    folderName = sprintf('fixedLatency_%sms', sanitize_threshold_text(stageB_settings.fixed_latency_threshold_ms));
end

function txt = sanitize_threshold_text(value)

    txt = regexprep(sprintf('%.3f', value), '0+$', '');
    txt = regexprep(txt, '\.$', '');
    txt = strrep(txt, '.', 'p');
end

function make_pvalue_distribution_plot(T, className, classOutDir, default_p_cutoff)

    validP = ~isnan(T.SALTPValue);
    p = T.SALTPValue(validP);
    if isempty(p)
        return
    end

    pPlot = max(p, 1e-12);
    pSorted = sort(pPlot);
    frac = (1:numel(pSorted))' / numel(pSorted);
    nSigDefault = nnz(p < default_p_cutoff);
    nSig005 = nnz(p < 0.05);
    fracSigDefault = nSigDefault / numel(p);
    fracSig005 = nSig005 / numel(p);

    fig = figure('Color', 'w');
    tiledlayout(3,1, 'Padding', 'compact', 'TileSpacing', 'compact');

    nexttile
    histogram(p, linspace(0, 1, 51), 'FaceColor', [0.5 0.5 0.5], 'EdgeColor', 'none');
    hold on
    xline(default_p_cutoff, 'r--', 'LineWidth', 1.5);
    xline(0.05, ':', 'Color', [0.4 0.4 0.4], 'LineWidth', 1.2);
    xline(0.001, ':', 'Color', [0.4 0.4 0.4], 'LineWidth', 1.2);
    xlabel('SALT-style p-value');
    ylabel('Cell count');
    title(sprintf('%s: SALT-style p-value distribution (full range)', className), 'Interpreter', 'none');
    grid on
    yLimits = ylim;
    xText = 0.62;
    yText = yLimits(2) * 0.92;
    text(xText, yText, sprintf('p < %.3g: %d / %d (%.1f%%)', default_p_cutoff, nSigDefault, numel(p), 100 * fracSigDefault), ...
        'FontSize', 9, 'BackgroundColor', 'w', 'Margin', 4);
    text(xText, yText - 0.10 * yLimits(2), sprintf('p < 0.05: %d / %d (%.1f%%)', nSig005, numel(p), 100 * fracSig005), ...
        'FontSize', 9, 'BackgroundColor', 'w', 'Margin', 4);

    nexttile
    histogram(p, linspace(0, 0.1, 51), 'FaceColor', [0.2 0.4 0.7], 'EdgeColor', 'none');
    hold on
    xline(default_p_cutoff, 'r--', 'LineWidth', 1.5);
    xline(0.05, ':', 'Color', [0.4 0.4 0.4], 'LineWidth', 1.2);
    xline(0.001, ':', 'Color', [0.4 0.4 0.4], 'LineWidth', 1.2);
    xlabel('SALT-style p-value');
    ylabel('Cell count');
    title(sprintf('%s: SALT-style p-value distribution (low range)', className), 'Interpreter', 'none');
    grid on

    nexttile
    semilogx(pSorted, frac, 'LineWidth', 2, 'Color', [0.2 0.4 0.7]);
    hold on
    xline(default_p_cutoff, 'r--', 'LineWidth', 1.5);
    xline(0.05, ':', 'Color', [0.4 0.4 0.4], 'LineWidth', 1.2);
    xline(0.001, ':', 'Color', [0.4 0.4 0.4], 'LineWidth', 1.2);
    xlabel('SALT-style p-value');
    ylabel('Cumulative fraction');
    title(sprintf('%s: cumulative SALT-style p-value distribution', className), 'Interpreter', 'none');
    grid on

    save_figure_png(fig, fullfile(classOutDir, sprintf('%s_SALT_PValueDistribution.png', className)));
    close(fig);
end

function make_tagged_latency_histogram(Tin, className, classOutDir, StageB, stageB_settings)

    stageB_settings = ensure_stage_b_settings_defaults(stageB_settings);
    taggedMask = Tin.StageA_Label == "tagged" & ~isnan(Tin.Latency_ms);
    if nnz(taggedMask) == 0
        return
    end

    latAll = Tin.Latency_ms(taggedMask);
    if is_integer_quantized(latAll)
        edges = (min(latAll) - 0.5):(max(latAll) + 0.5);
    else
        edges = compute_histogram_edges(latAll, stageB_settings.principal_latency_hist_bin_method);
    end
    lab = Tin.FinalLabel(taggedMask);

    figL = figure('Color', 'w');
    hold on

    histogram(latAll, edges, ...
        'FaceColor', [0.65 0.65 0.65], ...
        'FaceAlpha', 0.30, ...
        'EdgeColor', 'none', ...
        'DisplayName', 'all tagged');

    latDirect = latAll(lab == "direct");
    if ~isempty(latDirect)
        histogram(latDirect, edges, ...
            'FaceColor', [0.2 0.7 0.2], ...
            'FaceAlpha', 0.55, ...
            'EdgeColor', 'none', ...
            'DisplayName', 'direct');
    end

    latIndirect = latAll(lab == "indirect");
    if ~isempty(latIndirect)
        histogram(latIndirect, edges, ...
            'FaceColor', [0.85 0.33 0.10], ...
            'FaceAlpha', 0.65, ...
            'EdgeColor', 'none', ...
            'DisplayName', 'indirect');
    end

    if isfinite(StageB.LatencyThreshold_ms)
        xline(StageB.LatencyThreshold_ms, 'k--', 'LineWidth', 1.5, 'DisplayName', '7 ms cutoff');
    end

    xlabel('Latency (ms)');
    ylabel('Cell count');
    title(sprintf('%s: Tagged-cell latency histogram', className), 'Interpreter', 'none');
    legend('Location', 'best');
    grid on

    save_figure_png(figL, fullfile(classOutDir, sprintf('%s_Tagged_Latency_histogram.png', className)));
    close(figL);
end

function export_salt_manuscript_figure_outputs(T, StageB, outDir, exportSettings)

    if isempty(T) || height(T) == 0
        return
    end

    default_p_cutoff = get_export_setting_or_default(exportSettings, 'default_salt_cutoff', 0.01);
    min_reliability_for_tagged = get_export_setting_or_default(exportSettings, 'min_reliability_for_tagged', 0.10);
    analyze_dir = char(string(get_export_setting_or_default(exportSettings, 'analyze_dir', 'opto1')));
    delay_ms = get_export_setting_or_default(exportSettings, 'delay_ms', 0);
    response_window_ms = get_export_setting_or_default(exportSettings, 'response_window_ms', 15);
    stageB_settings = ensure_stage_b_settings_defaults(get_export_setting_or_default(exportSettings, 'stageB_settings', struct()));
    figureDir = fullfile(outDir, 'SALT-Figure');
    if ~exist(figureDir, 'dir')
        mkdir(figureDir);
    end

    n_showcase_examples_per_group = get_export_setting_or_default(exportSettings, 'n_showcase_examples_per_group', 3);
    indirect_example_latency_range_ms = get_export_setting_or_default(exportSettings, 'indirect_example_latency_range_ms', [9 14]);

    make_salt_workflow_schematic(figureDir, default_p_cutoff, min_reliability_for_tagged);
    make_manuscript_salt_pvalue_panel(T, figureDir, default_p_cutoff, 'all_cells');
    make_manuscript_tagged_latency_panel(T, figureDir, StageB, stageB_settings, 'all_cells');
    make_manuscript_final_count_panel(T, figureDir, 'all_cells');
    make_manuscript_composite_panel_set(T, figureDir, StageB, default_p_cutoff, stageB_settings);
    export_salt_example_panels(T, figureDir, analyze_dir, delay_ms, response_window_ms, ...
        min_reliability_for_tagged, n_showcase_examples_per_group, indirect_example_latency_range_ms);

    classNames = ["principal","interneuron"];
    for iClass = 1:numel(classNames)
        className = classNames(iClass);
        classMask = T.CellClass == className;
        if nnz(classMask) == 0
            continue
        end

        classTag = sanitize_filename(char(className));
        make_manuscript_salt_pvalue_panel(T(classMask,:), figureDir, default_p_cutoff, classTag);
        make_manuscript_tagged_latency_panel(T(classMask,:), figureDir, StageB, stageB_settings, classTag);
    end
end

function nameOut = sanitize_filename(nameIn)

    nameOut = regexprep(char(nameIn), '[^\w\s-]', '');
    nameOut = strtrim(nameOut);
    nameOut = regexprep(nameOut, '\s+', '_');
    if isempty(nameOut)
        nameOut = 'unnamed_graph';
    end
end

function d = kl_divergence(p, q)

    mask = p > 0 & q > 0;
    d = sum(p(mask) .* log2(p(mask) ./ q(mask)));
end

function StageB = initialize_salt_stage_b_struct(methodName)

    StageB = struct();
    StageB.Method = methodName;
    StageB.LatencyThreshold_ms = NaN;
    StageB.UsedLateLatencySplit = false;
    StageB.FixedThresholdSeparation = struct();
    StageB.RawMedians_Direct = struct('BaselineSD', NaN, 'Reliability', NaN, 'Latency_ms', NaN, 'Jitter_ms', NaN);
    StageB.RawMedians_Indirect = struct('BaselineSD', NaN, 'Reliability', NaN, 'Latency_ms', NaN, 'Jitter_ms', NaN);
    StageB.Diagnostics = struct();
    StageB.Diagnostics.NumTagged = NaN;
    StageB.Diagnostics.NumWithLatency = NaN;
    StageB.Diagnostics.ClassificationMode = "";
    StageB.Diagnostics.SharedAcrossClasses = true;
    StageB.Diagnostics.FixedLatencyThreshold_ms = NaN;
    StageB.Diagnostics.FixedThresholdOverrideApplied = false;
    StageB.Diagnostics.IndirectRule = "";
end

function StageB = populate_salt_stage_b_raw_medians(StageB, Tin, taggedMask)

    dirMask = taggedMask & (Tin.StageB_Label == "direct");
    indMask = taggedMask & (Tin.StageB_Label == "indirect");

    StageB.RawMedians_Direct.BaselineSD = median(Tin.BaselineSD(dirMask), 'omitnan');
    StageB.RawMedians_Direct.Reliability = median(Tin.Reliability(dirMask), 'omitnan');
    StageB.RawMedians_Direct.Latency_ms = median(Tin.Latency_ms(dirMask), 'omitnan');
    StageB.RawMedians_Direct.Jitter_ms = median(Tin.Jitter_ms(dirMask), 'omitnan');

    StageB.RawMedians_Indirect.BaselineSD = median(Tin.BaselineSD(indMask), 'omitnan');
    StageB.RawMedians_Indirect.Reliability = median(Tin.Reliability(indMask), 'omitnan');
    StageB.RawMedians_Indirect.Latency_ms = median(Tin.Latency_ms(indMask), 'omitnan');
    StageB.RawMedians_Indirect.Jitter_ms = median(Tin.Jitter_ms(indMask), 'omitnan');
end

function summaryStruct = summarize_fixed_threshold_separation(Tin, threshold_ms)

    if ~isfinite(threshold_ms)
        summaryStruct = struct( ...
            'Threshold_ms', NaN, ...
            'NumTaggedWithLatency', 0, ...
            'NumDirect', 0, ...
            'NumIndirect', 0, ...
            'MeanDirectLatency_ms', NaN, ...
            'MeanIndirectLatency_ms', NaN, ...
            'LatencyMeanGap_ms', NaN, ...
            'LatencyAshmanD', NaN, ...
            'MedianDirectLatency_ms', NaN, ...
            'MedianIndirectLatency_ms', NaN, ...
            'MedianDirectReliability', NaN, ...
            'MedianIndirectReliability', NaN, ...
            'MedianDirectJitter_ms', NaN, ...
            'MedianIndirectJitter_ms', NaN, ...
            'MedianDirectBaselineSD', NaN, ...
            'MedianIndirectBaselineSD', NaN, ...
            'PLatency', NaN, ...
            'ZLatency', NaN, ...
            'PReliability', NaN, ...
            'ZReliability', NaN, ...
            'PJitter', NaN, ...
            'ZJitter', NaN, ...
            'PBaselineSD', NaN, ...
            'ZBaselineSD', NaN);
        return
    end

    taggedMask = Tin.StageA_Label == "tagged" & ~isnan(Tin.Latency_ms);
    if nnz(taggedMask) == 0
        summaryStruct = summarize_fixed_threshold_separation(Tin, NaN);
        summaryStruct.Threshold_ms = threshold_ms;
        return
    end

    directMask = taggedMask & Tin.Latency_ms < threshold_ms;
    indirectMask = taggedMask & Tin.Latency_ms >= threshold_ms;

    [latency_p, latency_z] = run_ranksum_safe(Tin.Latency_ms(directMask), Tin.Latency_ms(indirectMask));
    [reliability_p, reliability_z] = run_ranksum_safe(Tin.Reliability(directMask), Tin.Reliability(indirectMask));
    [jitter_p, jitter_z] = run_ranksum_safe(Tin.Jitter_ms(directMask), Tin.Jitter_ms(indirectMask));
    [baseline_p, baseline_z] = run_ranksum_safe(Tin.BaselineSD(directMask), Tin.BaselineSD(indirectMask));

    meanDirectLatency = mean(Tin.Latency_ms(directMask), 'omitnan');
    meanIndirectLatency = mean(Tin.Latency_ms(indirectMask), 'omitnan');

    summaryStruct = struct( ...
        'Threshold_ms', threshold_ms, ...
        'NumTaggedWithLatency', nnz(taggedMask), ...
        'NumDirect', nnz(directMask), ...
        'NumIndirect', nnz(indirectMask), ...
        'MeanDirectLatency_ms', meanDirectLatency, ...
        'MeanIndirectLatency_ms', meanIndirectLatency, ...
        'LatencyMeanGap_ms', meanIndirectLatency - meanDirectLatency, ...
        'LatencyAshmanD', compute_ashman_d(Tin.Latency_ms(directMask), Tin.Latency_ms(indirectMask)), ...
        'MedianDirectLatency_ms', median(Tin.Latency_ms(directMask), 'omitnan'), ...
        'MedianIndirectLatency_ms', median(Tin.Latency_ms(indirectMask), 'omitnan'), ...
        'MedianDirectReliability', median(Tin.Reliability(directMask), 'omitnan'), ...
        'MedianIndirectReliability', median(Tin.Reliability(indirectMask), 'omitnan'), ...
        'MedianDirectJitter_ms', median(Tin.Jitter_ms(directMask), 'omitnan'), ...
        'MedianIndirectJitter_ms', median(Tin.Jitter_ms(indirectMask), 'omitnan'), ...
        'MedianDirectBaselineSD', median(Tin.BaselineSD(directMask), 'omitnan'), ...
        'MedianIndirectBaselineSD', median(Tin.BaselineSD(indirectMask), 'omitnan'), ...
        'PLatency', latency_p, ...
        'ZLatency', latency_z, ...
        'PReliability', reliability_p, ...
        'ZReliability', reliability_z, ...
        'PJitter', jitter_p, ...
        'ZJitter', jitter_z, ...
        'PBaselineSD', baseline_p, ...
        'ZBaselineSD', baseline_z);
end

function save_figure_png(figHandle, filename)
    [folder,~,~] = fileparts(filename);
    if ~exist(folder, 'dir')
        mkdir(folder);
    end
    exportgraphics(figHandle, filename, 'Resolution', 300);
end

function tf = is_integer_quantized(x)

    x = x(isfinite(x));
    if isempty(x)
        tf = false;
        return
    end

    tf = all(abs(x - round(x)) < 1e-9);
end

function edges = compute_histogram_edges(x, binMethod)

    x = x(isfinite(x));
    if isempty(x)
        edges = [0 1];
        return
    end

    if numel(unique(x)) == 1
        span = max(0.25, abs(x(1)) * 0.1 + eps);
        edges = [x(1) - span, x(1) + span];
        return
    end

    [~, edges] = histcounts(x, 'BinMethod', binMethod);

    if numel(edges) < 2 || any(~isfinite(edges)) || edges(1) == edges(end)
        edges = linspace(min(x), max(x), 16);
    end

    if numel(edges) < 2 || edges(1) == edges(end)
        span = max(std(x, 0, 'omitnan'), 0.25);
        edges = [min(x) - span, max(x) + span];
    end
end

function value = get_export_setting_or_default(export_settings, field_name, default_value)

    value = default_value;
    if isstruct(export_settings) && isfield(export_settings, field_name)
        candidate = export_settings.(field_name);
        if ~isempty(candidate)
            value = candidate;
        end
    end
end

function make_manuscript_salt_pvalue_panel(T, figureDir, default_p_cutoff, fileTag)

    validP = ~isnan(T.SALTPValue);
    pVals = T.SALTPValue(validP);
    if isempty(pVals)
        return
    end

    fig = figure('Color', 'w', 'Position', [100 100 520 420]);
    ax = axes(fig);
    plot_manuscript_salt_pvalue_histogram(ax, pVals, default_p_cutoff);

    xlabel(ax, 'SALT-style p-value');
    ylabel(ax, 'Cell count');
    title(ax, 'Step 1: Light-responsive cells', 'Interpreter', 'none');
    box(ax, 'off');

    yMax = ylim(ax);
    text(ax, default_p_cutoff, yMax(2) * 0.96, sprintf(' cutoff = %.3g', default_p_cutoff), ...
        'Color', [0.70 0.10 0.10], 'FontSize', 10, 'VerticalAlignment', 'top', ...
        'HorizontalAlignment', 'left', 'Interpreter', 'none');
    text(ax, max(xlim(ax)) * 0.98, yMax(2) * 0.96, sprintf('n = %d valid cells', numel(pVals)), ...
        'Color', [0.15 0.15 0.15], 'FontSize', 10, 'VerticalAlignment', 'top', ...
        'HorizontalAlignment', 'right', 'Interpreter', 'none');

    apply_manuscript_axes_style(ax);
    save_figure_outputs(fig, fullfile(figureDir, sprintf('Panel_Step1_SALT_PValueDistribution_%s', fileTag)));
    close(fig);
end

function make_manuscript_tagged_latency_panel(T, figureDir, StageB, stageB_settings, fileTag)

    taggedMask = T.StageA_Label == "tagged" & ~isnan(T.Latency_ms);
    if nnz(taggedMask) == 0
        return
    end

    latAll = T.Latency_ms(taggedMask);
    if is_integer_quantized(latAll)
        edges = (min(latAll) - 0.5):(max(latAll) + 0.5);
    else
        edges = compute_histogram_edges(latAll, stageB_settings.principal_latency_hist_bin_method);
    end

    labels = T.FinalLabel(taggedMask);
    latDirect = latAll(labels == "direct");
    latIndirect = latAll(labels == "indirect");

    fig = figure('Color', 'w', 'Position', [100 100 540 420]);
    ax = axes(fig);
    hold(ax, 'on');

    histogram(ax, latAll, edges, 'FaceColor', [0.75 0.75 0.75], 'FaceAlpha', 0.35, 'EdgeColor', 'none', 'DisplayName', 'all tagged');
    if ~isempty(latDirect)
        histogram(ax, latDirect, edges, 'FaceColor', [0.22 0.70 0.28], 'FaceAlpha', 0.60, 'EdgeColor', 'none', 'DisplayName', 'direct');
    end
    if ~isempty(latIndirect)
        histogram(ax, latIndirect, edges, 'FaceColor', [0.86 0.39 0.16], 'FaceAlpha', 0.70, 'EdgeColor', 'none', 'DisplayName', 'indirect');
    end

    if isfinite(StageB.LatencyThreshold_ms)
        xline(ax, StageB.LatencyThreshold_ms, 'k--', 'LineWidth', 1.5, 'DisplayName', '7 ms cutoff');
    end

    xlabel(ax, 'Latency (ms)');
    ylabel(ax, 'Cell count');
    title(ax, 'Step 2: Direct versus indirect tagged cells', 'Interpreter', 'none');
    legend(ax, 'Location', 'best');
    grid(ax, 'on');
    box(ax, 'off');

    yMax = ylim(ax);
    text(ax, max(xlim(ax)) * 0.98, yMax(2) * 0.96, sprintf('direct = %d\nindirect = %d', numel(latDirect), numel(latIndirect)), ...
        'Color', [0.15 0.15 0.15], 'FontSize', 10, 'VerticalAlignment', 'top', ...
        'HorizontalAlignment', 'right', 'Interpreter', 'none');

    apply_manuscript_axes_style(ax);
    save_figure_outputs(fig, fullfile(figureDir, sprintf('Panel_Step2_TaggedLatencyHistogram_%s', fileTag)));
    close(fig);
end

function make_manuscript_final_count_panel(T, figureDir, fileTag)

    classNames = ["principal","interneuron"];
    labelNames = ["untagged","indirect","direct"];
    countMatrix = zeros(numel(classNames) + 1, numel(labelNames));

    for iClass = 1:numel(classNames)
        classMask = T.CellClass == classNames(iClass);
        for iLabel = 1:numel(labelNames)
            countMatrix(iClass, iLabel) = nnz(classMask & T.FinalLabel == labelNames(iLabel));
        end
    end
    for iLabel = 1:numel(labelNames)
        countMatrix(end, iLabel) = nnz(T.FinalLabel == labelNames(iLabel));
    end

    fig = figure('Color', 'w', 'Position', [100 100 540 420]);
    ax = axes(fig);
    b = bar(ax, countMatrix, 'stacked', 'BarWidth', 0.68);
    b(1).FaceColor = [0.72 0.72 0.72];
    b(2).FaceColor = [0.86 0.39 0.16];
    b(3).FaceColor = [0.22 0.70 0.28];

    set(ax, 'XTick', 1:3, 'XTickLabel', {'principal', 'interneuron', 'all'});
    ylabel(ax, 'Cell count');
    title(ax, 'Two-step classification outcome', 'Interpreter', 'none');
    legend(ax, {'untagged', 'indirect', 'direct'}, 'Location', 'best');
    grid(ax, 'on');
    box(ax, 'off');

    apply_manuscript_axes_style(ax);
    save_figure_outputs(fig, fullfile(figureDir, sprintf('Panel_TwoStep_CountSummary_%s', fileTag)));
    close(fig);
end

function make_manuscript_composite_panel_set(T, figureDir, StageB, default_p_cutoff, stageB_settings)

    stageB_settings = ensure_stage_b_settings_defaults(stageB_settings);
    validP = ~isnan(T.SALTPValue);
    pVals = T.SALTPValue(validP);
    taggedMask = T.StageA_Label == "tagged" & ~isnan(T.Latency_ms);

    fig = figure('Color', 'w', 'Position', [100 100 1500 430]);
    tl = tiledlayout(fig, 1, 3, 'Padding', 'compact', 'TileSpacing', 'compact');

    ax1 = nexttile(tl, 1);
    if ~isempty(pVals)
        plot_manuscript_salt_pvalue_histogram(ax1, pVals, default_p_cutoff);
    end
    xlabel(ax1, 'SALT-style p-value');
    ylabel(ax1, 'Cell count');
    title(ax1, 'Step 1: Light-responsive cells', 'Interpreter', 'none');
    box(ax1, 'off');
    apply_manuscript_axes_style(ax1);

    ax2 = nexttile(tl, 2);
    if nnz(taggedMask) > 0
        latAll = T.Latency_ms(taggedMask);
        if is_integer_quantized(latAll)
            latencyEdges = (min(latAll) - 0.5):(max(latAll) + 0.5);
        else
            latencyEdges = compute_histogram_edges(latAll, stageB_settings.principal_latency_hist_bin_method);
        end
        labels = T.FinalLabel(taggedMask);
        histogram(ax2, latAll, latencyEdges, 'FaceColor', [0.75 0.75 0.75], 'FaceAlpha', 0.35, 'EdgeColor', 'none');
        hold(ax2, 'on');
        histogram(ax2, latAll(labels == "direct"), latencyEdges, 'FaceColor', [0.22 0.70 0.28], 'FaceAlpha', 0.60, 'EdgeColor', 'none');
        histogram(ax2, latAll(labels == "indirect"), latencyEdges, 'FaceColor', [0.86 0.39 0.16], 'FaceAlpha', 0.70, 'EdgeColor', 'none');
        if isfinite(StageB.LatencyThreshold_ms)
            xline(ax2, StageB.LatencyThreshold_ms, 'k--', 'LineWidth', 1.5);
        end
    end
    xlabel(ax2, 'Latency (ms)');
    ylabel(ax2, 'Cell count');
    title(ax2, 'Step 2: Direct versus indirect tagged cells', 'Interpreter', 'none');
    grid(ax2, 'on');
    box(ax2, 'off');
    apply_manuscript_axes_style(ax2);

    ax3 = nexttile(tl, 3);
    classNames = ["principal","interneuron"];
    labelNames = ["untagged","indirect","direct"];
    countMatrix = zeros(numel(classNames) + 1, numel(labelNames));
    for iClass = 1:numel(classNames)
        classMask = T.CellClass == classNames(iClass);
        for iLabel = 1:numel(labelNames)
            countMatrix(iClass, iLabel) = nnz(classMask & T.FinalLabel == labelNames(iLabel));
        end
    end
    for iLabel = 1:numel(labelNames)
        countMatrix(end, iLabel) = nnz(T.FinalLabel == labelNames(iLabel));
    end
    b = bar(ax3, countMatrix, 'stacked', 'BarWidth', 0.68);
    b(1).FaceColor = [0.72 0.72 0.72];
    b(2).FaceColor = [0.86 0.39 0.16];
    b(3).FaceColor = [0.22 0.70 0.28];
    set(ax3, 'XTick', 1:3, 'XTickLabel', {'principal', 'interneuron', 'all'});
    ylabel(ax3, 'Cell count');
    title(ax3, 'Two-step classification outcome', 'Interpreter', 'none');
    legend(ax3, {'untagged', 'indirect', 'direct'}, 'Location', 'best');
    grid(ax3, 'on');
    box(ax3, 'off');
    apply_manuscript_axes_style(ax3);

    save_figure_outputs(fig, fullfile(figureDir, 'SALT_MethodsFigure_PanelSet'));
    close(fig);
end

function make_salt_workflow_schematic(figureDir, default_p_cutoff, min_reliability_for_tagged)

    fig = figure('Color', 'w', 'Position', [100 100 1180 360]);
    ax = axes(fig, 'Position', [0 0 1 1], 'Visible', 'off');
    hold(ax, 'on');

    annotation(fig, 'textbox', [0.03 0.34 0.18 0.28], ...
        'String', sprintf('Optogenetic pulse train\n+\nspike timestamps'), ...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
        'FontName', 'Arial', 'FontSize', 12, 'LineWidth', 1.2, ...
        'BackgroundColor', [0.96 0.97 1.00], 'EdgeColor', [0.30 0.40 0.70]);

    annotation(fig, 'textbox', [0.28 0.20 0.22 0.56], ...
        'String', sprintf(['Step 1: SALT-style test\n' ...
        'First-spike latency distribution after laser onset\n' ...
        'versus baseline windows from the inter-pulse interval\n\n' ...
        'Tagged if:\n' ...
        'p < %.3g\n' ...
        'and reliability >= %.2f'], default_p_cutoff, min_reliability_for_tagged), ...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
        'FontName', 'Arial', 'FontSize', 12, 'LineWidth', 1.2, ...
        'BackgroundColor', [1.00 0.97 0.93], 'EdgeColor', [0.78 0.45 0.18]);

    annotation(fig, 'textbox', [0.58 0.20 0.22 0.56], ...
        'String', ['Step 2: fixed latency split' newline ...
        'Tagged cells only' newline ...
        'Direct if latency < 7 ms' newline ...
        'Indirect if latency >= 7 ms'], ...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
        'FontName', 'Arial', 'FontSize', 12, 'LineWidth', 1.2, ...
        'BackgroundColor', [0.94 0.99 0.94], 'EdgeColor', [0.22 0.60 0.22]);

    annotation(fig, 'textbox', [0.85 0.28 0.12 0.40], ...
        'String', ['Final labels' newline newline 'untagged' newline 'indirect' newline 'direct'], ...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
        'FontName', 'Arial', 'FontSize', 12, 'LineWidth', 1.2, ...
        'BackgroundColor', [0.97 0.97 0.97], 'EdgeColor', [0.45 0.45 0.45]);

    annotation(fig, 'arrow', [0.21 0.28], [0.48 0.48], 'LineWidth', 1.8);
    annotation(fig, 'arrow', [0.50 0.58], [0.48 0.48], 'LineWidth', 1.8);
    annotation(fig, 'arrow', [0.80 0.85], [0.48 0.48], 'LineWidth', 1.8);

    save_figure_outputs(fig, fullfile(figureDir, 'Panel_Workflow_Schematic'));
    close(fig);
end

function export_salt_example_panels(T, figureDir, analyze_dir, delay_ms, response_window_ms, ...
    min_reliability_for_tagged, nExamplesPerGroup, indirectLatencyRange_ms)

    sessionInfoPath = resolve_session_info_path();
    load(sessionInfoPath, 'sessInfo');

    exampleTable = select_salt_examples(T, min_reliability_for_tagged, nExamplesPerGroup, indirectLatencyRange_ms);
    if isempty(exampleTable) || height(exampleTable) == 0
        return
    end

    exampleDir = fullfile(figureDir, 'SelectedExamples');
    if ~exist(exampleDir, 'dir')
        mkdir(exampleDir);
    end

    for iExample = 1:height(exampleTable)
        exampleData = load_salt_example_data(exampleTable(iExample,:), sessInfo, analyze_dir, delay_ms, response_window_ms);
        if isempty(exampleData)
            continue
        end

        make_salt_example_figure(exampleData, fullfile(exampleDir, sprintf('Example_%s_%s_rank%d', ...
            sanitize_filename(char(exampleTable.CellClass(iExample))), ...
            sanitize_filename(char(exampleTable.FinalLabel(iExample))), ...
            exampleTable.ShowcaseRank(iExample))));
    end
end

function plot_manuscript_salt_pvalue_histogram(ax, pVals, default_p_cutoff)

    displayFloor = 1e-3;
    pVals = pVals(:);
    pValsDisplay = max(pVals, displayFloor);
    edges = sort(unique([logspace(log10(displayFloor), 0, 50), default_p_cutoff]));
    taggedMask = pVals < default_p_cutoff;

    hold(ax, 'on');
    histogram(ax, pValsDisplay(taggedMask), edges, ...
        'FaceColor', [0.22 0.46 0.74], 'EdgeColor', 'w', 'LineWidth', 0.5);
    histogram(ax, pValsDisplay(~taggedMask), edges, ...
        'FaceColor', [0.05 0.05 0.05], 'EdgeColor', 'w', 'LineWidth', 0.5);

    set(ax, 'XScale', 'log');
    xline(ax, default_p_cutoff, 'r--', 'LineWidth', 1.6);
    xlim(ax, [displayFloor 1]);
    xticks(ax, [displayFloor default_p_cutoff 1]);
    xticklabels(ax, {'0', '0.01', '1'});
    ax.XGrid = 'off';
    ax.YGrid = 'off';
    ax.XMinorGrid = 'off';
    ax.YMinorGrid = 'off';
end

function save_figure_outputs(figHandle, fileBase, panelRows, applyManuscriptFormatting)
    if nargin < 3 || isempty(panelRows)
        panelRows = 1;
    end
    if nargin < 4 || isempty(applyManuscriptFormatting)
        applyManuscriptFormatting = true;
    end

    if applyManuscriptFormatting
        prepare_manuscript_figure_for_export(figHandle, panelRows);
    end
    save_figure_png(figHandle, sprintf('%s.png', fileBase));
    save_figure_pdf(figHandle, sprintf('%s.pdf', fileBase));
end

function apply_manuscript_axes_style(ax)
    set(ax, 'FontName', 'Arial', 'FontSize', 14, 'FontWeight', 'bold', 'LineWidth', 1.0, 'TickDir', 'out');
end

function exampleTable = select_salt_examples(T, min_reliability_for_tagged, nExamplesPerGroup, indirectLatencyRange_ms)

    if nargin < 3 || isempty(nExamplesPerGroup)
        nExamplesPerGroup = 3;
    end
    if nargin < 4 || isempty(indirectLatencyRange_ms) || numel(indirectLatencyRange_ms) < 2
        indirectLatencyRange_ms = [9 14];
    end

    rows = zeros(0,1);
    showcaseRank = zeros(0,1);
    selectionClass = strings(0,1);
    selectionLabel = strings(0,1);

    [manualRows, manualMessages] = resolve_manual_direct_examples(T);
    for iMessage = 1:numel(manualMessages)
        fprintf('%s\n', manualMessages{iMessage});
    end
    if ~isempty(manualRows)
        rows = [rows; manualRows(:)]; %#ok<AGROW>
        showcaseRank = [showcaseRank; (1:numel(manualRows))']; %#ok<AGROW>
        selectionClass = [selectionClass; T.CellClass(manualRows)]; %#ok<AGROW>
        selectionLabel = [selectionLabel; repmat("direct_manual", numel(manualRows), 1)]; %#ok<AGROW>
    end

    classOrder = ["principal","interneuron"];
    for iClass = 1:numel(classOrder)
        className = classOrder(iClass);

        indirectMask = T.CellClass == className & T.FinalLabel == "indirect" & ...
            ~isnan(T.SALTPValue) & ~isnan(T.Reliability) & isfinite(T.Latency_ms);
        indirectMask = indirectMask & T.Reliability >= min_reliability_for_tagged;
        indirectIdx = find(indirectMask);
        if isempty(indirectIdx)
            indirectIdx = find(T.CellClass == className & T.FinalLabel == "indirect" & ...
                ~isnan(T.SALTPValue) & ~isnan(T.Reliability) & isfinite(T.Latency_ms));
        end

        if ~isempty(indirectIdx)
            indirectLatency = T.Latency_ms(indirectIdx);
            preferredMask = indirectLatency >= indirectLatencyRange_ms(1) & indirectLatency <= indirectLatencyRange_ms(2);
            if any(preferredMask)
                indirectIdx = indirectIdx(preferredMask);
                indirectLatency = indirectLatency(preferredMask);
            end

            latencyCenter = mean(indirectLatencyRange_ms);
            distanceFromCenter = abs(indirectLatency - latencyCenter);
            score = [-log10(max(T.SALTPValue(indirectIdx), 1e-12)), T.Reliability(indirectIdx), -distanceFromCenter, indirectLatency];
            [~, order] = sortrows(score, [-1 -2 -3 -4]);

            nTake = min(nExamplesPerGroup, numel(order));
            chosenRows = indirectIdx(order(1:nTake));
            rows = [rows; chosenRows(:)]; %#ok<AGROW>
            showcaseRank = [showcaseRank; (1:nTake)']; %#ok<AGROW>
            selectionClass = [selectionClass; repmat(className, nTake, 1)]; %#ok<AGROW>
            selectionLabel = [selectionLabel; repmat("indirect", nTake, 1)]; %#ok<AGROW>
        end

        labelName = "untagged";
        mask = T.CellClass == className & T.FinalLabel == labelName & ...
            ~isnan(T.SALTPValue) & ~isnan(T.Reliability);
        mask = mask & T.Reliability >= min_reliability_for_tagged;

        idx = find(mask);
        if isempty(idx)
            idx = find(T.CellClass == className & T.FinalLabel == labelName & ...
                ~isnan(T.SALTPValue) & ~isnan(T.Reliability));
        end

        if isempty(idx)
            continue
        end

        score = [T.SALTPValue(idx), T.Reliability(idx)];
        [~, order] = sortrows(score, [-1 -2]);

        nTake = min(nExamplesPerGroup, numel(order));
        chosenRows = idx(order(1:nTake));
        rows = [rows; chosenRows(:)]; %#ok<AGROW>
        showcaseRank = [showcaseRank; (1:nTake)']; %#ok<AGROW>
        selectionClass = [selectionClass; repmat(className, nTake, 1)]; %#ok<AGROW>
        selectionLabel = [selectionLabel; repmat(labelName, nTake, 1)]; %#ok<AGROW>
    end

    if isempty(rows)
        exampleColumns = salt_example_export_columns(T);
        exampleTable = T([], exampleColumns);
        exampleTable.ShowcaseRank = zeros(0,1);
        exampleTable.SelectionClass = strings(0,1);
        exampleTable.SelectionLabel = strings(0,1);
        return
    end

    exampleColumns = salt_example_export_columns(T);
    exampleTable = T(rows, exampleColumns);
    exampleTable.SelectionClass = selectionClass;
    exampleTable.SelectionLabel = selectionLabel;
    exampleTable.ShowcaseRank = showcaseRank;
end

function columns = salt_example_export_columns(T)

    columns = {'SessionIndex','INumber','Cell','CellClass','FinalLabel','SALTPValue', ...
        'Reliability','Latency_ms','Jitter_ms','RoomID','AnimalID'};
    if ismember('RoomName', T.Properties.VariableNames)
        columns = {'SessionIndex','INumber','Cell','CellClass','FinalLabel','SALTPValue', ...
            'Reliability','Latency_ms','Jitter_ms','RoomName','RoomID','AnimalID'};
    end
end

function exampleData = load_salt_example_data(exampleRow, sessInfo, analyze_dir, delay_ms, response_window_ms)

    exampleData = struct([]);
    sessionIdx = exampleRow.SessionIndex;
    cellIdx = exampleRow.Cell;

    if sessionIdx > numel(sessInfo) || ~isfield(sessInfo(sessionIdx), 'mainDir') || ~isfield(sessInfo(sessionIdx), 'tList')
        return
    end

    mainDir = normalize_text_field(sessInfo(sessionIdx).mainDir);
    tListName = normalize_text_field(sessInfo(sessionIdx).tList);
    if isempty(mainDir) || isempty(tListName)
        return
    end

    optoDir = fullfile(mainDir, analyze_dir);
    saltMetricPath = fullfile(optoDir, sprintf('SALTMetric_%dmsDelay.mat', delay_ms));
    if ~exist(saltMetricPath, 'file')
        return
    end

    S = load(saltMetricPath, 'sessionMetrics');
    if ~isfield(S, 'sessionMetrics')
        return
    end
    sessionMetrics = S.sessionMetrics;

    ttlOn = [];
    ttlOff = [];
    baselineWindows = [];
    if isfield(sessionMetrics, 'ttl_on')
        ttlOn = sessionMetrics.ttl_on(:);
    end
    if isfield(sessionMetrics, 'ttl_off')
        ttlOff = sessionMetrics.ttl_off(:);
    end
    if isfield(sessionMetrics, 'baseline_windows')
        baselineWindows = sessionMetrics.baseline_windows;
    end

    if isempty(ttlOn) || isempty(ttlOff)
        return
    end

    if isfield(sessionMetrics, 'tt_files') && numel(sessionMetrics.tt_files) >= cellIdx
        ttFile = sessionMetrics.tt_files{cellIdx};
    else
        tt_files = read_tt_list(fullfile(mainDir, tListName));
        if numel(tt_files) < cellIdx
            return
        end
        ttFile = tt_files{cellIdx};
    end

    try
        spikeData = readSpikeDataOnly(optoDir, {ttFile});
        spikeTimes = fixSpikes(spikeData);
    catch
        return
    end

    if isempty(spikeTimes) || isempty(spikeTimes{1})
        cellSpikes = [];
    else
        cellSpikes = spikeTimes{1}(:);
    end

    pulseWidth_ms = median((ttlOff - ttlOn) * 1000, 'omitnan');
    if ~isfinite(pulseWidth_ms) || pulseWidth_ms <= 0
        pulseWidth_ms = 10;
    end

    exampleData = struct( ...
        'SessionIndex', sessionIdx, ...
        'INumber', exampleRow.INumber, ...
        'Cell', cellIdx, ...
        'CellClass', exampleRow.CellClass, ...
        'FinalLabel', exampleRow.FinalLabel, ...
        'SALTPValue', exampleRow.SALTPValue, ...
        'Reliability', exampleRow.Reliability, ...
        'Latency_ms', exampleRow.Latency_ms, ...
        'Jitter_ms', exampleRow.Jitter_ms, ...
        'AnimalID', exampleRow.AnimalID, ...
        'RoomID', exampleRow.RoomID, ...
        'CellSpikes', cellSpikes, ...
        'TTLon', ttlOn, ...
        'BaselineWindows', baselineWindows, ...
        'ResponseWindow_ms', response_window_ms, ...
        'PulseWidth_ms', pulseWidth_ms);
end

function make_salt_example_figure(exampleData, fileBase)

    display_window_ms = 25;
    pre_ms = display_window_ms;
    post_ms = display_window_ms;
    stim_window_ms = exampleData.ResponseWindow_ms;
    example_panel_width_scale = 0.60;
    example_panel_width_px = round(520 * example_panel_width_scale);
    example_panel_height_px = 700;

    [stimRasterTimes_ms, stimRasterTrials] = collect_aligned_spikes_ms(exampleData.CellSpikes, exampleData.TTLon, pre_ms / 1000, post_ms / 1000);
    [bin_centers_ms, counts_all_per_sec, counts_first_per_sec] = compute_example_focused_psth( ...
        exampleData.CellSpikes, exampleData.TTLon, exampleData.ResponseWindow_ms / 1000, exampleData.ResponseWindow_ms / 1000);
    nLatencyBins = max(1, round(exampleData.ResponseWindow_ms));
    stimCats = first_spike_categories(exampleData.CellSpikes, exampleData.TTLon, exampleData.ResponseWindow_ms / 1000, 0.001, nLatencyBins);
    stimProb = category_histogram(stimCats, nLatencyBins + 1);
    if ~isempty(exampleData.BaselineWindows)
        baseCats = first_spike_categories(exampleData.CellSpikes, exampleData.BaselineWindows(:,1), exampleData.ResponseWindow_ms / 1000, 0.001, nLatencyBins);
        baseProb = category_histogram(baseCats, nLatencyBins + 1);
    else
        baseProb = nan(1, nLatencyBins + 1);
    end

    fig = figure('Color', 'w', 'Position', [100 100 example_panel_width_px example_panel_height_px]);
    tl = tiledlayout(fig, 3, 1, 'Padding', 'compact', 'TileSpacing', 'compact');

    ax1 = nexttile(tl, 1);
    if ~isempty(stimRasterTimes_ms)
        scatter(ax1, stimRasterTimes_ms, stimRasterTrials, 10, [0.1 0.1 0.1], 'filled');
    end
    hold(ax1, 'on');
    patch(ax1, [0 stim_window_ms stim_window_ms 0], [0 0 max(1, numel(exampleData.TTLon) + 1) max(1, numel(exampleData.TTLon) + 1)], ...
        [0.75 0.88 1.00], 'EdgeColor', 'none', 'FaceAlpha', 0.4);
    if ~isempty(stimRasterTimes_ms)
        scatter(ax1, stimRasterTimes_ms, stimRasterTrials, 10, [0.1 0.1 0.1], 'filled');
    end
    set(ax1, 'YDir', 'reverse', 'XLim', [-pre_ms post_ms]);
    ylabel(ax1, 'Pulse #');
    title(ax1, sprintf('%s %s example | animal=%s | i=%d cell=%d | p=%s', ...
        char(exampleData.CellClass), char(exampleData.FinalLabel), char(string(exampleData.AnimalID)), exampleData.INumber, exampleData.Cell, format_p_value(exampleData.SALTPValue)), ...
        'Interpreter', 'none');
    box(ax1, 'off');
    apply_example_axes_style(ax1);

    ax2 = nexttile(tl, 2);
    hold(ax2, 'on');
    y_max = max([double(counts_all_per_sec(:)); double(counts_first_per_sec(:))]);
    if ~isfinite(y_max) || y_max <= 0
        y_max = 1;
    end
    patch(ax2, [0 stim_window_ms stim_window_ms 0], [0 0 y_max * 1.08 y_max * 1.08], ...
        [0.75 0.88 1.00], 'EdgeColor', 'none', 'FaceAlpha', 0.35);
    bar(ax2, bin_centers_ms, double(counts_all_per_sec(:)), 0.9, ...
        'FaceColor', [0.72 0.72 0.72], 'EdgeColor', 'none');
    bar(ax2, bin_centers_ms, double(counts_first_per_sec(:)), 0.58, ...
        'FaceColor', [0.05 0.05 0.05], 'EdgeColor', 'none');
    plot(ax2, [0 0], [0 y_max * 1.08], '--', 'Color', [0.2 0.2 0.2], 'LineWidth', 1);
    plot(ax2, [stim_window_ms stim_window_ms], [0 y_max * 1.08], '--', 'Color', [0.2 0.2 0.2], 'LineWidth', 1);
    xlim(ax2, [-display_window_ms display_window_ms]);
    ylim(ax2, [0 y_max * 1.08]);
    xlabel(ax2, 'Time from pulse onset (ms)');
    ylabel(ax2, 'Firing rate (Hz)');
    title(ax2, 'Opto1 PSTH (-25 to +25 ms, rounded 1 ms bins)', 'Interpreter', 'none');
    legend(ax2, {'Stim window', 'All spikes (grey)', 'First spikes (black)'}, 'Location', 'northwest', 'Box', 'off');
    text(ax2, bin_centers_ms(end) * 0.98, y_max * 0.96, sprintf('Rel = %.2f\nLat = %.2f ms', exampleData.Reliability, exampleData.Latency_ms), ...
        'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', 'Interpreter', 'none');
    box(ax2, 'off');
    apply_example_axes_style(ax2);

    ax3 = nexttile(tl, 3);
    xVals = 1:(nLatencyBins + 1);
    xLabels = [arrayfun(@(k) sprintf('%d-%d', k - 1, k), 1:nLatencyBins, 'UniformOutput', false), {'no spike'}];
    if all(isfinite(baseProb))
        plot(ax3, xVals, baseProb, '-o', 'Color', [0.55 0.55 0.55], 'LineWidth', 1.4, ...
            'MarkerFaceColor', [0.70 0.70 0.70], 'DisplayName', 'baseline');
        hold(ax3, 'on');
    else
        hold(ax3, 'on');
    end
    plot(ax3, xVals, stimProb, '-o', 'Color', [0.82 0.30 0.10], 'LineWidth', 1.6, ...
        'MarkerFaceColor', [0.86 0.39 0.16], 'DisplayName', 'laser');
    xline(ax3, nLatencyBins + 0.5, ':', 'Color', [0.55 0.55 0.55], 'LineWidth', 1.0, 'HandleVisibility', 'off');
    set(ax3, 'XTick', xVals, 'XTickLabel', xLabels, 'XTickLabelRotation', 45);
    xlim(ax3, [1 nLatencyBins + 1]);
    xlabel(ax3, 'First-spike category (ms)');
    ylabel(ax3, 'Fraction of windows');
    title(ax3, 'SALT first-spike category distribution', 'Interpreter', 'none');
    legend(ax3, 'Location', 'best');
    box(ax3, 'off');
    apply_example_axes_style(ax3);

    save_figure_outputs(fig, fileBase, [], false);
    close(fig);
end

function [pVal, zVal] = run_ranksum_safe(group1, group2)

    pVal = NaN;
    zVal = NaN;

    group1 = group1(isfinite(group1));
    group2 = group2(isfinite(group2));

    if isempty(group1) || isempty(group2)
        return
    end

    try
        [pVal, ~, stats] = ranksum(group1, group2);
        if isfield(stats, 'zval')
            zVal = stats.zval;
        end
    catch
        pVal = NaN;
        zVal = NaN;
    end
end

function ashmanD = compute_ashman_d(group1, group2)

    group1 = group1(isfinite(group1));
    group2 = group2(isfinite(group2));

    if isempty(group1) || isempty(group2)
        ashmanD = NaN;
        return
    end

    mu1 = mean(group1, 'omitnan');
    mu2 = mean(group2, 'omitnan');
    sd1 = std(group1, 0, 'omitnan');
    sd2 = std(group2, 0, 'omitnan');
    denom = sqrt(sd1^2 + sd2^2);

    if ~isfinite(denom) || denom <= 0
        ashmanD = NaN;
        return
    end

    ashmanD = sqrt(2) * abs(mu2 - mu1) / denom;
end

function save_figure_pdf(figHandle, filename)
    [folder,~,~] = fileparts(filename);
    if ~exist(folder, 'dir')
        mkdir(folder);
    end
    exportgraphics(figHandle, filename, 'ContentType', 'vector');
end

function prepare_manuscript_figure_for_export(figHandle, panelRows)

    if nargin < 2 || isempty(panelRows)
        panelRows = 1;
    end

    targetPanelHeightCm = 5;
    targetFigureHeightCm = max(panelRows, 1) * targetPanelHeightCm;

    originalUnits = get(figHandle, 'Units');
    set(figHandle, 'Units', 'pixels');
    posPx = get(figHandle, 'Position');

    if numel(posPx) < 4 || ~isfinite(posPx(3)) || ~isfinite(posPx(4)) || posPx(3) <= 0 || posPx(4) <= 0
        aspectRatio = 1;
    else
        aspectRatio = posPx(3) / posPx(4);
    end

    set(figHandle, 'Units', 'centimeters');
    posCm = get(figHandle, 'Position');
    posCm(3) = max(targetFigureHeightCm * aspectRatio, 2.5);
    posCm(4) = targetFigureHeightCm;
    set(figHandle, 'Position', posCm, 'PaperPositionMode', 'auto');
    set(figHandle, 'Units', originalUnits);

    apply_manuscript_figure_font_style(figHandle);
    drawnow;
end

function [rows, messages] = resolve_manual_direct_examples(T)

specs = manual_direct_example_specs();
rows = zeros(0, 1);
messages = cell(0, 1);

for iSpec = 1:numel(specs)
    spec = specs(iSpec);
    matchMask = T.INumber == spec.INumber & T.Cell == spec.Cell;

    if strlength(spec.AnimalID) > 0
        matchMask = matchMask & T.AnimalID == spec.AnimalID;
    end

    idx = find(matchMask);
    if isempty(idx)
        messages{end + 1, 1} = sprintf('[SALT examples] Manual direct example not found: animal=%s i=%d cell=%d %s', ...
            char(spec.AnimalID), spec.INumber, spec.Cell, char(spec.TTFile)); %#ok<AGROW>
        continue
    end

    directIdx = idx(T.FinalLabel(idx) == "direct");
    if ~isempty(directIdx)
        idx = directIdx;
    end

    rows(end + 1, 1) = idx(1); %#ok<AGROW>
end

end

function [bin_centers_ms, counts_all_per_sec, counts_first_per_sec] = compute_example_focused_psth(cell_spikes, ttl_on, first_spike_window_s, stim_duration_s)

bin_sec = 0.001;
window_sec = 0.025;
bin_centers_ms = (-25:25)';
counts_all = zeros(size(bin_centers_ms));
counts_first = zeros(size(bin_centers_ms));

if isempty(cell_spikes) || isempty(ttl_on)
    counts_all_per_sec = counts_all;
    counts_first_per_sec = counts_first;
    return
end

if ~isfinite(first_spike_window_s) || first_spike_window_s <= 0
    first_spike_window_s = 0.015;
end
if ~isfinite(stim_duration_s) || stim_duration_s <= 0
    stim_duration_s = first_spike_window_s;
end

for pulse_idx = 1:numel(ttl_on)
    ttl_time = ttl_on(pulse_idx);
    keep_mask = cell_spikes >= (ttl_time - window_sec) & cell_spikes <= (ttl_time + window_sec);
    relative_spikes = cell_spikes(keep_mask) - ttl_time;
    if isempty(relative_spikes)
        continue
    end

    rounded_all_ms = round(relative_spikes .* 1000);
    valid_all = rounded_all_ms >= bin_centers_ms(1) & rounded_all_ms <= bin_centers_ms(end);
    rounded_all_ms = rounded_all_ms(valid_all);
    if ~isempty(rounded_all_ms)
        counts_all = counts_all + accumarray(rounded_all_ms(:) - bin_centers_ms(1) + 1, 1, [numel(bin_centers_ms), 1]);
    end

    first_candidates = relative_spikes(relative_spikes >= 0 & relative_spikes <= first_spike_window_s);
    if isempty(first_candidates)
        continue
    end

    first_latency_ms = round(first_candidates(1) * 1000);
    if first_latency_ms >= 0 && first_latency_ms <= round(first_spike_window_s * 1000)
        counts_first(first_latency_ms - bin_centers_ms(1) + 1) = counts_first(first_latency_ms - bin_centers_ms(1) + 1) + 1;
    end
end

counts_all_per_sec = counts_all ./ numel(ttl_on) ./ bin_sec;
counts_first_per_sec = counts_first ./ numel(ttl_on) ./ bin_sec;

end

function [times_ms, trials] = collect_aligned_spikes_ms(spikeTimes, eventTimes, pre_s, post_s)

    times_ms = [];
    trials = [];
    if isempty(spikeTimes) || isempty(eventTimes)
        return
    end

    for iEvent = 1:numel(eventTimes)
        relTimes = spikeTimes - eventTimes(iEvent);
        keepMask = relTimes >= -pre_s & relTimes <= post_s;
        if any(keepMask)
            rel_kept = relTimes(keepMask) * 1000;
            times_ms = [times_ms; rel_kept(:)]; %#ok<AGROW>
            trials = [trials; repmat(iEvent, numel(rel_kept), 1)]; %#ok<AGROW>
        end
    end
end

function txt = format_p_value(pVal)

    if ~isfinite(pVal)
        txt = 'n/a';
    elseif pVal < 1e-4
        txt = '< 1e-4';
    else
        txt = sprintf('%.3g', pVal);
    end
end

function apply_example_axes_style(ax)
    set(ax, 'FontName', 'Arial', 'FontSize', 11, 'FontWeight', 'normal', 'LineWidth', 1.0, 'TickDir', 'out');
end

function apply_manuscript_figure_font_style(figHandle)

    fontNamedObjects = findall(figHandle, '-property', 'FontName');
    for iObj = 1:numel(fontNamedObjects)
        set(fontNamedObjects(iObj), 'FontName', 'Arial');
    end

    fontSizedObjects = findall(figHandle, '-property', 'FontSize');
    for iObj = 1:numel(fontSizedObjects)
        set(fontSizedObjects(iObj), 'FontSize', 14);
    end

    fontWeightedObjects = findall(figHandle, '-property', 'FontWeight');
    for iObj = 1:numel(fontWeightedObjects)
        set(fontWeightedObjects(iObj), 'FontWeight', 'bold');
    end
end

function specs = manual_direct_example_specs()

specs = [ ...
    struct('AnimalID', "1715", 'INumber', 150, 'Cell', 1,  'TTFile', "TT1_01"); ...
    struct('AnimalID', "1715", 'INumber', 147, 'Cell', 5,  'TTFile', "TT1_05"); ...
    struct('AnimalID', "1715", 'INumber', 143, 'Cell', 4,  'TTFile', "TT2_01"); ...
    struct('AnimalID', "1713", 'INumber', 253, 'Cell', 21, 'TTFile', "TT4_01"); ...
    struct('AnimalID', "524",  'INumber', 32,  'Cell', 16, 'TTFile', "TT3_01"); ...
    struct('AnimalID', "521",  'INumber', 23,  'Cell', 2,  'TTFile', "TT1_02"); ...
    struct('AnimalID', "521",  'INumber', 24,  'Cell', 1,  'TTFile', "TT1_01"); ...
    struct('AnimalID', "1715", 'INumber', 153, 'Cell', 10, 'TTFile', "TT3_02")];

end
