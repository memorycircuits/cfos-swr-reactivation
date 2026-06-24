function CSV_file_export_first(custom_settings)
% CSV_file_export_first
% Export curated CSV files and non-SALT comparison outputs from cached SALT classification.
%
% Run SALT_opto_analysis first to create SALT_TwoStep_ClassificationResults.mat.

close all
clc
rng(1)

if nargin < 1 || isempty(custom_settings)
    custom_settings = struct();
end

preferredOutDir = get_override_value(custom_settings, 'preferredOutDir', fullfile(pwd, 'OptoMetricComparison_SALT'));
preferredCsvExportDir = get_override_value(custom_settings, 'preferredCsvExportDir', fullfile(pwd, 'CSVfiles', 'curatedFiles'));
export_only_animals_with_direct_cells = get_override_value(custom_settings, 'export_only_animals_with_direct_cells', true);
stageB_settings = ensure_stage_b_settings_defaults(struct('fixed_latency_threshold_ms', 7));

outDir = resolve_output_dir(preferredOutDir, stageB_settings);
csvExportDir = resolve_csv_export_dir(preferredCsvExportDir, stageB_settings);
allCellsPath = char(string(get_override_value(custom_settings, 'allCellsPath', '')));
if isempty(allCellsPath)
    allCellsPath = resolve_all_cells_path();
end
resultsCachePath = char(string(get_override_value(custom_settings, 'resultsCachePath', '')));
if isempty(resultsCachePath)
    resultsCachePath = fullfile(outDir, 'SALT_TwoStep_ClassificationResults.mat');
end

fprintf('\n=== CSV Export From SALT Classification ===\n');
fprintf('SALT cache: %s\n', resultsCachePath);
fprintf('Output directory: %s\n', outDir);
fprintf('Curated CSV directory: %s\n', csvExportDir);

if ~exist(resultsCachePath, 'file')
    error('Could not find cached SALT classification results at %s. Run SALT_opto_analysis first.', resultsCachePath);
end

loadedAllCells = load(allCellsPath, 'All_Cells_combined');
All_Cells_combined = loadedAllCells.All_Cells_combined;
cachedResults = load(resultsCachePath, 'T', 'SaltTwoStep', 'ExportSettings');
if ~isfield(cachedResults, 'T') || ~isfield(cachedResults, 'SaltTwoStep') || ~isfield(cachedResults, 'ExportSettings')
    error('Cached SALT classification results are incomplete: %s', resultsCachePath);
end

cachedResults.ExportSettings.export_only_animals_with_direct_cells = export_only_animals_with_direct_cells;
cachedResults.ExportSettings.stageB_settings = stageB_settings;

run_curated_exports(cachedResults.T, cachedResults.SaltTwoStep, ...
    cachedResults.ExportSettings, outDir, csvExportDir, All_Cells_combined);
fprintf('\nSaved curated CSV exports to:\n%s\n', csvExportDir);

end

function run_curated_exports(T, SaltTwoStep, ExportSettings, outDir, csvExportDir, All_Cells_combined)

    if isempty(T) || height(T) == 0
        fprintf('No valid cells available for SALT comparison.\n');
        return
    end

    T = add_place_cell_flag_to_table(T, All_Cells_combined);

    default_salt_cutoff = ExportSettings.default_salt_cutoff;
    min_reliability_for_tagged = get_export_setting_or_default(ExportSettings, 'min_reliability_for_tagged', 0.10);
    old_sd_cutoff = get_export_setting_or_default(ExportSettings, 'old_sd_cutoff', 5);
    export_only_animals_with_direct_cells = get_export_setting_or_default(ExportSettings, 'export_only_animals_with_direct_cells', true);
    stageB_settings = ensure_stage_b_settings_defaults(ExportSettings.stageB_settings);

    if export_only_animals_with_direct_cells
        [T, includedAnimalIds] = restrict_exports_to_direct_animals(T);
        fprintf('Curated SALT exports restricted to %d animals with at least one direct cell.\n', numel(includedAnimalIds));
        if isempty(T)
            fprintf('No animals met the direct-cell inclusion rule for curated exports.\n');
            return
        end
    else
        fprintf('Curated SALT exports include all animals in the SALT table.\n');
    end

    write_main_room_context_place_field_counts_excel(T, outDir);

    classes_to_run = ["principal","interneuron","place cells","non-place principal"];
    include_in_combined_results = [true,true,false,false];
    Results_all = table;
    Summary = table;
    FeatureSummary = table;

    for ic = 1:numel(classes_to_run)
        className = classes_to_run(ic);
        mask = resolve_class_export_mask(T, className);

        if nnz(mask) == 0
            continue
        end

        classOutDir = fullfile(outDir, sanitize_filename(char(className)));
        if ~exist(classOutDir, 'dir')
            mkdir(classOutDir);
        end

        [Results_class, Thresholds] = run_old_activation_metric( ...
            T(mask,:), char(className), old_sd_cutoff);

        Results_class.NewActivated_Default = compute_salt_activation_mask( ...
            Results_class, default_salt_cutoff, min_reliability_for_tagged);
        validInputs = get_salt_activation_input_mask(Results_class);

        Results_class.CompareValid = Results_class.OldMetricValid & validInputs;
        Results_class.ComparisonLabel = repmat("missing", height(Results_class), 1);

        bothMask = Results_class.CompareValid & Results_class.OldActivated & Results_class.NewActivated_Default;
        oldOnlyMask = Results_class.CompareValid & Results_class.OldActivated & ~Results_class.NewActivated_Default;
        newOnlyMask = Results_class.CompareValid & ~Results_class.OldActivated & Results_class.NewActivated_Default;
        neitherMask = Results_class.CompareValid & ~Results_class.OldActivated & ~Results_class.NewActivated_Default;

        Results_class.ComparisonLabel(neitherMask) = "neither";
        Results_class.ComparisonLabel(oldOnlyMask) = "old only";
        Results_class.ComparisonLabel(newOnlyMask) = "new only";
        Results_class.ComparisonLabel(bothMask) = "both";

        SummaryRow = summarize_class_comparison(Results_class, className, default_salt_cutoff, min_reliability_for_tagged, Thresholds);
        if include_in_combined_results(ic)
            Results_all = [Results_all; Results_class]; %#ok<AGROW>
        end
        Summary = [Summary; SummaryRow]; %#ok<AGROW>

        writetable(Results_class, fullfile(classOutDir, sprintf('%s_SALTMetricComparison.csv', className)));
        write_class_summary_txt(SummaryRow, fullfile(classOutDir, sprintf('%s_SALTMetricComparisonSummary.txt', className)));

        roomGroups = ["main","control"];
        for iRoom = 1:numel(roomGroups)
            roomGroupName = roomGroups(iRoom);
            roomMask = Results_class.RoomGroup == roomGroupName;
            if nnz(roomMask) == 0
                continue
            end

            roomOutDir = fullfile(classOutDir, sprintf('%s_sessions', roomGroupName));
            if ~exist(roomOutDir, 'dir')
                mkdir(roomOutDir);
            end

            FeatureSummary_room = compare_direct_vs_rest_feature_sets( ...
                Results_class(roomMask,:), char(className), char(roomGroupName), roomOutDir);
            if ~isempty(FeatureSummary_room)
                FeatureSummary = [FeatureSummary; FeatureSummary_room]; %#ok<AGROW>
            end
        end

        labelGroups = ["direct","non_direct"];
        for iLabel = 1:numel(labelGroups)
            labelName = labelGroups(iLabel);
            FeatureSummary_label = compare_room_groups_within_label_feature_sets( ...
                Results_class, char(className), char(labelName), classOutDir);
            if ~isempty(FeatureSummary_label)
                FeatureSummary = [FeatureSummary; FeatureSummary_label]; %#ok<AGROW>
            end

            FeatureSummary_specificRooms = compare_main_vs_specific_control_rooms_swr_features( ...
                Results_class, char(className), char(labelName), classOutDir);
            if ~isempty(FeatureSummary_specificRooms)
                FeatureSummary = [FeatureSummary; FeatureSummary_specificRooms]; %#ok<AGROW>
            end
        end

        FeatureSummary_fourCondition = compare_main_control_four_condition_feature_sets( ...
            Results_class, char(className), classOutDir, csvExportDir);
        if ~isempty(FeatureSummary_fourCondition)
            FeatureSummary = [FeatureSummary; FeatureSummary_fourCondition]; %#ok<AGROW>
        end

        FeatureSummary_fourConditionSession = compare_main_control_four_condition_aggregated_feature_sets( ...
            Results_class, char(className), classOutDir, 'session');
        if ~isempty(FeatureSummary_fourConditionSession)
            FeatureSummary = [FeatureSummary; FeatureSummary_fourConditionSession]; %#ok<AGROW>
        end

        FeatureSummary_fourConditionAnimal = compare_main_control_four_condition_aggregated_feature_sets( ...
            Results_class, char(className), classOutDir, 'animal');
        if ~isempty(FeatureSummary_fourConditionAnimal)
            FeatureSummary = [FeatureSummary; FeatureSummary_fourConditionAnimal]; %#ok<AGROW>
        end
    end

    if ~isempty(Results_all)
        allOutDir = fullfile(outDir, 'all_classes');
        if ~exist(allOutDir, 'dir')
            mkdir(allOutDir);
        end

        Results_all.NewActivated_Default = compute_salt_activation_mask( ...
            Results_all, default_salt_cutoff, min_reliability_for_tagged);
        validInputs_all = get_salt_activation_input_mask(Results_all);
        Results_all.CompareValid = Results_all.OldMetricValid & validInputs_all;
        Results_all.ComparisonLabel = repmat("missing", height(Results_all), 1);
        Results_all.ComparisonLabel(Results_all.CompareValid & ~Results_all.OldActivated & ~Results_all.NewActivated_Default) = "neither";
        Results_all.ComparisonLabel(Results_all.CompareValid & Results_all.OldActivated & ~Results_all.NewActivated_Default) = "old only";
        Results_all.ComparisonLabel(Results_all.CompareValid & ~Results_all.OldActivated & Results_all.NewActivated_Default) = "new only";
        Results_all.ComparisonLabel(Results_all.CompareValid & Results_all.OldActivated & Results_all.NewActivated_Default) = "both";

        CombinedThresholds = struct('ScoreThreshold', old_sd_cutoff);
        CombinedSummary = summarize_class_comparison(Results_all, "all_classes", default_salt_cutoff, min_reliability_for_tagged, CombinedThresholds);
        writetable(Results_all, fullfile(outDir, 'AllClasses_SALTMetricComparison.csv'));
        writetable(Summary, fullfile(outDir, 'AllClasses_SALTMetricComparisonSummary.csv'));
        write_class_summary_txt(CombinedSummary, fullfile(allOutDir, 'all_classes_SALTMetricComparisonSummary.txt'));
    end

    if ~isempty(FeatureSummary)
        writetable(FeatureSummary, fullfile(outDir, 'SALT_DirectVsRest_SWRFeatureSummary.csv'));
    end

    if all(ismember({'ThetaPhasePreferred','ThetaPhaseP','ThetaPhaseNSpikes'}, T.Properties.VariableNames))
        try
            theta_preferred_phase_export_dotplots( ...
                'T', T, ...
                'All_Cells_combined', All_Cells_combined, ...
                'csvExportDir', csvExportDir);
        catch ME
            warning('Theta preferred phase export failed: %s', ME.message);
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

function outDir = resolve_csv_export_dir(preferredOutDir, stageB_settings)

    if nargin < 1 || isempty(preferredOutDir)
        preferredOutDir = fullfile(pwd, 'CSVfiles', 'curatedFiles');
    end

    outDir = preferredOutDir;
    parentDir = fileparts(preferredOutDir);

    if ~(isempty(parentDir) || exist(parentDir, 'dir'))
        outDir = fullfile(pwd, 'CSVfiles', 'curatedFiles');
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

function [T_out, includedAnimalIds] = restrict_exports_to_direct_animals(T_in)

    includedAnimalIds = unique(T_in.AnimalID(T_in.FinalLabel == "direct" & T_in.AnimalID ~= ""), 'stable');
    if isempty(includedAnimalIds)
        T_out = T_in([],:);
        return
    end

    includeMask = ismember(T_in.AnimalID, includedAnimalIds);
    T_out = T_in(includeMask, :);
end

function T = add_place_cell_flag_to_table(T, All_Cells_combined)

    T.IsPlaceCell = false(height(T), 1);
    numericColumns = { ...
        'SpeedScore', 'SpatialInfo_A', ...
        'AverageRate_OF1', 'AverageRate_OF2', 'AverageRate_OF3', 'AverageRate_OFMean', ...
        'TemporalPeakRate', 'TemporalPeakRate_OF1', 'TemporalPeakRate_OF2', 'TemporalPeakRate_OF3', ...
        'SpatialPeakRate', 'SpatialPeakRate_OF1', 'SpatialPeakRate_OF2', 'SpatialPeakRate_OF3', 'SpatialPeakRate_OFMax', ...
        'PlaceFieldNumber', 'PlaceFieldNumber_OF1', 'PlaceFieldNumber_OF2', 'PlaceFieldNumber_OF3'};
    T = ensure_numeric_table_columns(T, numericColumns);

    for s = 1:numel(All_Cells_combined)
        sessionMask = T.SessionIndex == s;
        if ~any(sessionMask)
            continue
        end

        placeMask = [];
        spatialClassNumeric = get_field_or_empty(All_Cells_combined(s), 'final_classification_numeric');
        spatialClassText = get_field_or_empty(All_Cells_combined(s), 'final_cells_classification');
        speedScore = get_field_or_empty(All_Cells_combined(s), 'of_avg_speedScore');
        spatialInfo_A = get_field_or_empty(All_Cells_combined(s), 'context_A_spatial_info');
        peakRate = get_first_available_field(All_Cells_combined(s), {'peak_rate', 'of_avg_peak_fr_rate'});
        avgRateOF = get_first_available_field(All_Cells_combined(s), {'of_avg_fir_rate', 'classific_firingRate'});
        avgRate1 = get_first_available_field(All_Cells_combined(s), {'averRate_of1', 'AverageRate_OF1', 'meanRate_of1', 'of1_fir_rate', 'fir_rate_of1'});
        avgRate2 = get_first_available_field(All_Cells_combined(s), {'averRate_of2', 'AverageRate_OF2', 'meanRate_of2', 'of2_fir_rate', 'fir_rate_of2'});
        avgRate3 = get_first_available_field(All_Cells_combined(s), {'averRate_of3', 'AverageRate_OF3', 'meanRate_of3', 'of3_fir_rate', 'fir_rate_of3'});
        temporalPeakRate1 = get_temporal_peak_rate_of(All_Cells_combined(s), 1);
        temporalPeakRate2 = get_temporal_peak_rate_of(All_Cells_combined(s), 2);
        temporalPeakRate3 = get_temporal_peak_rate_of(All_Cells_combined(s), 3);
        spatialPeak1 = get_first_available_field(All_Cells_combined(s), {'spatialPeakRate_of1', 'SpatialPeakRate_OF1', 'peakRate_of1'});
        spatialPeak2 = get_first_available_field(All_Cells_combined(s), {'spatialPeakRate_of2', 'SpatialPeakRate_OF2', 'peakRate_of2'});
        spatialPeak3 = get_first_available_field(All_Cells_combined(s), {'spatialPeakRate_of3', 'SpatialPeakRate_OF3', 'peakRate_of3'});
        spatialPeakFallback = get_first_available_field(All_Cells_combined(s), {'spatialPeakRate', 'SpatialPeakRate', 'spatial_peak_rate'});
        pfnum1 = get_first_available_field(All_Cells_combined(s), {'PF_fieldNumbers_of1', 'of1_place_field_numbers'});
        pfnum2 = get_first_available_field(All_Cells_combined(s), {'PF_fieldNumbers_of2', 'of2_place_field_numbers'});
        pfnum3 = get_first_available_field(All_Cells_combined(s), {'PF_fieldNumbers_of3', 'of3_place_field_numbers'});

        if ~isempty(spatialClassNumeric)
            placeMask = double(spatialClassNumeric(:)) == 2;
        elseif ~isempty(spatialClassText)
            spatialClassText = string(spatialClassText(:));
            placeMask = contains(lower(strtrim(spatialClassText)), "place");
        end

        sessionRows = find(sessionMask);
        sessionRows = sessionRows(:);
        cellIdx = T.Cell(sessionMask);
        cellIdx = cellIdx(:);

        if ~isempty(speedScore)
            validSpeedIdx = isfinite(cellIdx) & cellIdx >= 1 & cellIdx <= numel(speedScore);
            T.SpeedScore(sessionRows(validSpeedIdx)) = speedScore(cellIdx(validSpeedIdx));
        end
        if ~isempty(spatialInfo_A)
            validSpatialInfoIdx = isfinite(cellIdx) & cellIdx >= 1 & cellIdx <= numel(spatialInfo_A);
            T.SpatialInfo_A(sessionRows(validSpatialInfoIdx)) = spatialInfo_A(cellIdx(validSpatialInfoIdx));
        end
        validCellIdx = cellIdx(isfinite(cellIdx) & cellIdx > 0);
        maxCellIdx = max([0; validCellIdx(:)]);
        nMetric = max([maxCellIdx, numel(peakRate), numel(avgRateOF), ...
            numel(avgRate1), numel(avgRate2), numel(avgRate3), ...
            numel(temporalPeakRate1), numel(temporalPeakRate2), numel(temporalPeakRate3), ...
            numel(spatialPeak1), numel(spatialPeak2), numel(spatialPeak3), numel(spatialPeakFallback), ...
            numel(pfnum1), numel(pfnum2), numel(pfnum3), 0]);

        if nMetric > 0
            peak_rate_pad = pad_numeric_vector(peakRate, nMetric);
            avg_rate_of_pad = pad_numeric_vector(avgRateOF, nMetric);
            avg_rate1_pad = pad_numeric_vector(avgRate1, nMetric);
            avg_rate2_pad = pad_numeric_vector(avgRate2, nMetric);
            avg_rate3_pad = pad_numeric_vector(avgRate3, nMetric);
            temporal_peak_rate1_pad = pad_numeric_vector(temporalPeakRate1, nMetric);
            temporal_peak_rate2_pad = pad_numeric_vector(temporalPeakRate2, nMetric);
            temporal_peak_rate3_pad = pad_numeric_vector(temporalPeakRate3, nMetric);
            spatial_peak1_pad = pad_numeric_vector(spatialPeak1, nMetric);
            spatial_peak2_pad = pad_numeric_vector(spatialPeak2, nMetric);
            spatial_peak3_pad = pad_numeric_vector(spatialPeak3, nMetric);
            spatial_peak_fallback_pad = pad_numeric_vector(spatialPeakFallback, nMetric);
            pfnum1_pad = pad_numeric_vector(pfnum1, nMetric);
            pfnum2_pad = pad_numeric_vector(pfnum2, nMetric);
            pfnum3_pad = pad_numeric_vector(pfnum3, nMetric);

            avg_rate1_export = row_nanmean_with_fallback(avg_rate1_pad, avg_rate_of_pad);
            avg_rate2_export = row_nanmean_with_fallback(avg_rate2_pad, avg_rate_of_pad);
            avg_rate3_export = row_nanmean_with_fallback(avg_rate3_pad, avg_rate_of_pad);
            avg_rate_of_mean = row_nanmean_with_fallback([avg_rate1_pad avg_rate2_pad avg_rate3_pad], avg_rate_of_pad);
            temporal_peak_rate1_export = row_nanmean_with_fallback(temporal_peak_rate1_pad, peak_rate_pad);
            temporal_peak_rate2_export = row_nanmean_with_fallback(temporal_peak_rate2_pad, peak_rate_pad);
            temporal_peak_rate3_export = row_nanmean_with_fallback(temporal_peak_rate3_pad, peak_rate_pad);
            spatial_peak_overall = row_nanmax_with_fallback( ...
                [spatial_peak1_pad spatial_peak2_pad spatial_peak3_pad], spatial_peak_fallback_pad);
            pfnum_mean = row_nanmean([pfnum1_pad pfnum2_pad pfnum3_pad]);

            T = assign_metric_by_cell_index(T, 'TemporalPeakRate', sessionRows, cellIdx, peak_rate_pad);
            T = assign_metric_by_cell_index(T, 'AverageRate_OF1', sessionRows, cellIdx, avg_rate1_export);
            T = assign_metric_by_cell_index(T, 'AverageRate_OF2', sessionRows, cellIdx, avg_rate2_export);
            T = assign_metric_by_cell_index(T, 'AverageRate_OF3', sessionRows, cellIdx, avg_rate3_export);
            T = assign_metric_by_cell_index(T, 'AverageRate_OFMean', sessionRows, cellIdx, avg_rate_of_mean);
            T = assign_metric_by_cell_index(T, 'TemporalPeakRate_OF1', sessionRows, cellIdx, temporal_peak_rate1_export);
            T = assign_metric_by_cell_index(T, 'TemporalPeakRate_OF2', sessionRows, cellIdx, temporal_peak_rate2_export);
            T = assign_metric_by_cell_index(T, 'TemporalPeakRate_OF3', sessionRows, cellIdx, temporal_peak_rate3_export);
            T = assign_metric_by_cell_index(T, 'SpatialPeakRate', sessionRows, cellIdx, spatial_peak_overall);
            T = assign_metric_by_cell_index(T, 'SpatialPeakRate_OF1', sessionRows, cellIdx, spatial_peak1_pad);
            T = assign_metric_by_cell_index(T, 'SpatialPeakRate_OF2', sessionRows, cellIdx, spatial_peak2_pad);
            T = assign_metric_by_cell_index(T, 'SpatialPeakRate_OF3', sessionRows, cellIdx, spatial_peak3_pad);
            T = assign_metric_by_cell_index(T, 'SpatialPeakRate_OFMax', sessionRows, cellIdx, spatial_peak_overall);
            T = assign_metric_by_cell_index(T, 'PlaceFieldNumber', sessionRows, cellIdx, pfnum_mean);
            T = assign_metric_by_cell_index(T, 'PlaceFieldNumber_OF1', sessionRows, cellIdx, pfnum1_pad);
            T = assign_metric_by_cell_index(T, 'PlaceFieldNumber_OF2', sessionRows, cellIdx, pfnum2_pad);
            T = assign_metric_by_cell_index(T, 'PlaceFieldNumber_OF3', sessionRows, cellIdx, pfnum3_pad);
        end

        if isempty(placeMask)
            continue
        end

        validIdx = isfinite(cellIdx) & cellIdx >= 1 & cellIdx <= numel(placeMask);
        T.IsPlaceCell(sessionRows(validIdx)) = placeMask(cellIdx(validIdx));
    end
end

function [Tin, Thresholds] = run_old_activation_metric(Tin, className, old_sd_cutoff)

    n = height(Tin);
    Tin.OldMetricValue = Tin.BaselineSD;
    Tin.OldActivated = false(n, 1);
    Tin.OldMetricValid = ~isnan(Tin.BaselineSD);

    Thresholds = struct();
    Thresholds.CellClass = className;
    Thresholds.ScoreThreshold = old_sd_cutoff;
    Thresholds.MetricName = 'SD_0_ms_delay';

    validIdx = find(Tin.OldMetricValid);
    Tin.OldActivated(validIdx) = Tin.BaselineSD(validIdx) > old_sd_cutoff;
end

function SummaryRow = summarize_class_comparison(T, className, default_p_cutoff, min_reliability_for_tagged, Thresholds)

    compareMask = T.CompareValid;
    oldPos = T.OldActivated(compareMask);
    newPos = T.NewActivated_Default(compareMask);

    both = nnz(oldPos & newPos);
    oldOnly = nnz(oldPos & ~newPos);
    newOnly = nnz(~oldPos & newPos);
    neither = nnz(~oldPos & ~newPos);

    unionCount = both + oldOnly + newOnly;
    if unionCount == 0
        jaccard = NaN;
    else
        jaccard = both / unionCount;
    end

    if nnz(newPos) == 0
        precision = NaN;
    else
        precision = both / nnz(newPos);
    end

    if nnz(oldPos) == 0
        recall = NaN;
    else
        recall = both / nnz(oldPos);
    end

    agreement = mean(oldPos == newPos, 'omitnan');

    SummaryRow = table( ...
        string(className), ...
        height(T), ...
        nnz(T.OldMetricValid), ...
        nnz(~isnan(T.SALTPValue)), ...
        nnz(compareMask), ...
        nnz(T.OldActivated & compareMask), ...
        nnz(T.NewActivated_Default & compareMask), ...
        both, oldOnly, newOnly, neither, ...
        jaccard, precision, recall, agreement, ...
        Thresholds.ScoreThreshold, ...
        default_p_cutoff, ...
        min_reliability_for_tagged, ...
        'VariableNames', {'CellClass','TotalCells','OldMetricValid','PValueValid','CompareValid', ...
        'OldActivated','NewActivated','Both','OldOnly','NewOnly','Neither', ...
        'Jaccard','Precision','Recall','Agreement','OldMetricThreshold','PValueCutoff','ReliabilityCutoff'});
end

function write_class_summary_txt(SummaryRow, filename)

    fid = fopen(filename, 'w');
    if fid == -1
        return
    end
    cleaner = onCleanup(@() fclose(fid)); %#ok<NASGU>

    fprintf(fid, 'Cell class: %s\n', char(string(SummaryRow.CellClass(1))));
    fprintf(fid, 'Total cells: %d\n', SummaryRow.TotalCells);
    fprintf(fid, 'Cells with old metric: %d\n', SummaryRow.OldMetricValid);
    fprintf(fid, 'Cells with SALT p-value: %d\n', SummaryRow.PValueValid);
    fprintf(fid, 'Cells compared: %d\n', SummaryRow.CompareValid);
    fprintf(fid, 'Old activated: %d\n', SummaryRow.OldActivated);
    fprintf(fid, 'SALT activated (default cutoff): %d\n', SummaryRow.NewActivated);
    fprintf(fid, 'Both: %d\n', SummaryRow.Both);
    fprintf(fid, 'Old only: %d\n', SummaryRow.OldOnly);
    fprintf(fid, 'SALT only: %d\n', SummaryRow.NewOnly);
    fprintf(fid, 'Neither: %d\n', SummaryRow.Neither);
    fprintf(fid, 'Jaccard: %.4f\n', SummaryRow.Jaccard);
    fprintf(fid, 'Precision: %.4f\n', SummaryRow.Precision);
    fprintf(fid, 'Recall: %.4f\n', SummaryRow.Recall);
    fprintf(fid, 'Agreement: %.4f\n', SummaryRow.Agreement);
    fprintf(fid, 'Old metric threshold (SD_0_ms_delay): %.4f\n', SummaryRow.OldMetricThreshold);
    fprintf(fid, 'Default SALT cutoff: %.4f\n', SummaryRow.PValueCutoff);
    fprintf(fid, 'Reliability cutoff for SALT tagging: %.4f\n', SummaryRow.ReliabilityCutoff);
end

function validMask = get_salt_activation_input_mask(T)

    validMask = ~isnan(T.SALTPValue) & ~isnan(T.Reliability);
end

function activationMask = compute_salt_activation_mask(T, p_cutoff, min_reliability_for_tagged)

    validMask = get_salt_activation_input_mask(T);
    activationMask = false(height(T), 1);
    activationMask(validMask) = T.SALTPValue(validMask) < p_cutoff & ...
        T.Reliability(validMask) >= min_reliability_for_tagged;
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

function SummaryTable = compare_direct_vs_rest_feature_sets(T, className, roomGroupName, classOutDir)

    [swrDefs, generalDefs] = build_feature_metric_defs(className);
    swrDefs = filter_metric_defs_by_available_data(T, swrDefs);
    generalDefs = filter_metric_defs_by_available_data(T, generalDefs);
    SummaryTable = table;

    figSWR = create_metric_tiled_figure(size(swrDefs, 1));
    for i = 1:size(swrDefs, 1)
        nexttile
        statsRow = plot_direct_vs_rest_metric_panel(gca, T, swrDefs{i,1}, swrDefs{i,2}, swrDefs{i,3}, i == 1);
        SummaryTable = [SummaryTable; statsRow]; %#ok<AGROW>
    end
    sgtitle(sprintf('%s %s sessions: untagged + indirect vs direct (SALT)', className, roomGroupName), 'Interpreter', 'none');
    save_figure_png(figSWR, fullfile(classOutDir, sprintf('%s_%s_SALT_DirectVsRest_SWRFeatures.png', className, roomGroupName)));
    close(figSWR);

    figGeneral = create_metric_tiled_figure(size(generalDefs, 1));
    for i = 1:size(generalDefs, 1)
        nexttile
        statsRow = plot_direct_vs_rest_metric_panel(gca, T, generalDefs{i,1}, generalDefs{i,2}, generalDefs{i,3}, i == 1);
        SummaryTable = [SummaryTable; statsRow]; %#ok<AGROW>
    end
    sgtitle(sprintf('%s %s sessions: untagged + indirect vs direct (SALT)', className, roomGroupName), 'Interpreter', 'none');
    save_figure_png(figGeneral, fullfile(classOutDir, sprintf('%s_%s_SALT_DirectVsRest_OtherFeatures.png', className, roomGroupName)));
    close(figGeneral);

    if ~isempty(SummaryTable)
        SummaryTable.CellClass = repmat(string(className), height(SummaryTable), 1);
        SummaryTable.ComparisonType = repmat("direct_vs_rest", height(SummaryTable), 1);
        SummaryTable.Subset = repmat(string(roomGroupName), height(SummaryTable), 1);
        SummaryTable.ControlRoom = repmat("", height(SummaryTable), 1);
        SummaryTable = movevars(SummaryTable, {'CellClass','ComparisonType','Subset','ControlRoom'}, 'Before', 1);
        writetable(SummaryTable, fullfile(classOutDir, sprintf('%s_%s_SALT_DirectVsRest_FeatureSummary.csv', className, roomGroupName)));
    end
end

function SummaryTable = compare_room_groups_within_label_feature_sets(T, className, labelName, classOutDir)

    [swrDefs, generalDefs] = build_feature_metric_defs(className);
    swrDefs = filter_metric_defs_by_available_data(T, swrDefs);
    generalDefs = filter_metric_defs_by_available_data(T, generalDefs);
    SummaryTable = table;
    labelMask = resolve_label_mask(T, labelName) & (T.RoomGroup == "main" | T.RoomGroup == "control");
    labelDisplayName = label_display_name(labelName);

    if nnz(labelMask) == 0
        return
    end

    roomOutDir = fullfile(classOutDir, 'room_comparisons');
    if ~exist(roomOutDir, 'dir')
        mkdir(roomOutDir);
    end

    figSWR = create_metric_tiled_figure(size(swrDefs, 1));
    for i = 1:size(swrDefs, 1)
        nexttile
        statsRow = plot_room_group_metric_panel(gca, T, swrDefs{i,1}, swrDefs{i,2}, swrDefs{i,3}, labelName);
        SummaryTable = [SummaryTable; statsRow]; %#ok<AGROW>
    end
    sgtitle(sprintf('%s %s cells: main-room vs control rooms', className, labelDisplayName), 'Interpreter', 'none');
    save_figure_png(figSWR, fullfile(roomOutDir, sprintf('%s_%s_SALT_MainVsControl_SWRFeatures.png', className, labelName)));
    close(figSWR);

    figGeneral = create_metric_tiled_figure(size(generalDefs, 1));
    for i = 1:size(generalDefs, 1)
        nexttile
        statsRow = plot_room_group_metric_panel(gca, T, generalDefs{i,1}, generalDefs{i,2}, generalDefs{i,3}, labelName);
        SummaryTable = [SummaryTable; statsRow]; %#ok<AGROW>
    end
    sgtitle(sprintf('%s %s cells: main-room vs control rooms', className, labelDisplayName), 'Interpreter', 'none');
    save_figure_png(figGeneral, fullfile(roomOutDir, sprintf('%s_%s_SALT_MainVsControl_OtherFeatures.png', className, labelName)));
    close(figGeneral);

    if ~isempty(SummaryTable)
        SummaryTable.CellClass = repmat(string(className), height(SummaryTable), 1);
        SummaryTable.ComparisonType = repmat("main_vs_control", height(SummaryTable), 1);
        SummaryTable.Subset = repmat(string(labelName), height(SummaryTable), 1);
        SummaryTable.ControlRoom = repmat("All Control Rooms", height(SummaryTable), 1);
        SummaryTable = movevars(SummaryTable, {'CellClass','ComparisonType','Subset','ControlRoom'}, 'Before', 1);
        writetable(SummaryTable, fullfile(roomOutDir, sprintf('%s_%s_SALT_MainVsControl_FeatureSummary.csv', className, labelName)));
    end
end

function SummaryTable = compare_main_vs_specific_control_rooms_swr_features(T, className, labelName, classOutDir)

    [swrDefs, ~] = build_feature_metric_defs(className);
    swrDefs = filter_metric_defs_by_available_data(T, swrDefs);
    SummaryTable = table;
    labelDisplayName = label_display_name(labelName);
    roomDefs = { ...
        19, 'control_room_familiar', 'Control Room Familiar'; ...
        20, 'control_room_novel', 'Control Room Novel'};
    roomOutDir = fullfile(classOutDir, 'room_comparisons');

    if ~exist(roomOutDir, 'dir')
        mkdir(roomOutDir);
    end

    for iRoom = 1:size(roomDefs, 1)
        controlRoomId = roomDefs{iRoom, 1};
        controlRoomFileName = roomDefs{iRoom, 2};
        controlRoomDisplayName = roomDefs{iRoom, 3};

        labelMask = resolve_label_mask(T, labelName);
        compareMask = labelMask & (T.RoomID == 15 | T.RoomID == controlRoomId);
        if nnz(compareMask) == 0
            continue
        end

        figSWR = create_metric_tiled_figure(size(swrDefs, 1));
        SummaryTable_room = table;

        for i = 1:size(swrDefs, 1)
            nexttile
            statsRow = plot_specific_room_metric_panel(gca, T, swrDefs{i,1}, swrDefs{i,2}, swrDefs{i,3}, labelName, controlRoomId, controlRoomDisplayName);
            SummaryTable_room = [SummaryTable_room; statsRow]; %#ok<AGROW>
        end

        sgtitle(sprintf('%s %s cells: main-room vs %s', className, labelDisplayName, controlRoomDisplayName), 'Interpreter', 'none');
        save_figure_png(figSWR, fullfile(roomOutDir, sprintf('%s_%s_SALT_MainVs%s_SWRFeatures.png', className, labelName, controlRoomFileName)));
        close(figSWR);

        if ~isempty(SummaryTable_room)
            SummaryTable_room.CellClass = repmat(string(className), height(SummaryTable_room), 1);
            SummaryTable_room.ComparisonType = repmat("main_vs_specific_control", height(SummaryTable_room), 1);
            SummaryTable_room.Subset = repmat(string(labelName), height(SummaryTable_room), 1);
            SummaryTable_room.ControlRoom = repmat(string(controlRoomDisplayName), height(SummaryTable_room), 1);
            SummaryTable_room = movevars(SummaryTable_room, {'CellClass','ComparisonType','Subset','ControlRoom'}, 'Before', 1);
            writetable(SummaryTable_room, fullfile(roomOutDir, sprintf('%s_%s_SALT_MainVs%s_SWRFeatureSummary.csv', className, labelName, controlRoomFileName)));
            SummaryTable = [SummaryTable; SummaryTable_room]; %#ok<AGROW>
        end
    end
end

function SummaryTable = compare_main_control_four_condition_feature_sets(T, className, classOutDir, csvExportDir)

    [swrDefs, generalDefs] = build_feature_metric_defs(className);
    swrDefs = filter_metric_defs_by_available_data(T, swrDefs);
    generalDefs = filter_metric_defs_by_available_data(T, generalDefs);
    SummaryTable = table;
    roomOutDir = fullfile(classOutDir, 'room_comparisons');

    if ~exist(roomOutDir, 'dir')
        mkdir(roomOutDir);
    end

    figSWR = create_metric_tiled_figure(size(swrDefs, 1));
    for i = 1:size(swrDefs, 1)
        nexttile
        statsRows = plot_four_condition_metric_panel(gca, T, swrDefs{i,1}, swrDefs{i,2}, swrDefs{i,3}, i == 1);
        SummaryTable = [SummaryTable; statsRows]; %#ok<AGROW>
        export_four_condition_csv_columns(T, swrDefs{i,1}, swrDefs{i,2}, swrDefs{i,3}, className, csvExportDir);
    end
    sgtitle(sprintf('%s cells: main-room/control cfos vs non-fos', className), 'Interpreter', 'none');
    save_figure_png(figSWR, fullfile(roomOutDir, sprintf('%s_SALT_MainControl_FourCondition_SWRFeatures.png', className)));
    close(figSWR);

    figGeneral = create_metric_tiled_figure(size(generalDefs, 1));
    for i = 1:size(generalDefs, 1)
        nexttile
        statsRows = plot_four_condition_metric_panel(gca, T, generalDefs{i,1}, generalDefs{i,2}, generalDefs{i,3}, i == 1);
        SummaryTable = [SummaryTable; statsRows]; %#ok<AGROW>
        export_four_condition_csv_columns(T, generalDefs{i,1}, generalDefs{i,2}, generalDefs{i,3}, className, csvExportDir);
    end
    sgtitle(sprintf('%s cells: main-room/control cfos vs non-fos', className), 'Interpreter', 'none');
    save_figure_png(figGeneral, fullfile(roomOutDir, sprintf('%s_SALT_MainControl_FourCondition_OtherFeatures.png', className)));
    close(figGeneral);

    if ~isempty(SummaryTable)
        SummaryTable.CellClass = repmat(string(className), height(SummaryTable), 1);
        SummaryTable.ComparisonType = repmat("main_control_four_condition", height(SummaryTable), 1);
        SummaryTable.Subset = repmat("cfos_vs_non_fos", height(SummaryTable), 1);
        SummaryTable.ControlRoom = repmat("All Control Rooms", height(SummaryTable), 1);
        SummaryTable = movevars(SummaryTable, {'CellClass','ComparisonType','Subset','ControlRoom'}, 'Before', 1);
        writetable(SummaryTable, fullfile(roomOutDir, sprintf('%s_SALT_MainControl_FourCondition_FeatureSummary.csv', className)));
    end
end

function write_main_room_context_place_field_counts_excel(T, outDir)

    requiredColumns = {'SessionIndex', 'INumber', 'Cell', 'AnimalID', ...
        'RoomID', 'CellClass', 'FinalLabel', ...
        'PlaceFieldNumber_OF1', 'PlaceFieldNumber_OF2', 'PlaceFieldNumber_OF3'};
    if ~all(ismember(requiredColumns, T.Properties.VariableNames))
        warning('main-room place-field count export skipped because required columns are missing.');
        return
    end

    if ~exist(outDir, 'dir')
        mkdir(outDir);
    end

    [IncludedValues, contextNames, contextLabels] = build_main_room_context_place_field_values(T);
    PooledCounts = summarize_main_room_context_place_field_values(IncludedValues, contextNames, contextLabels, false);
    ByContextCounts = summarize_main_room_context_place_field_values(IncludedValues, contextNames, contextLabels, true);
    Notes = main_room_context_place_field_notes_table();

    outputFile = fullfile(outDir, 'main_room_context_place_field_counts.xlsx');

    try
        if exist(outputFile, 'file')
            delete(outputFile);
        end
        writetable(Notes, outputFile, 'Sheet', 'Notes');
        writetable(PooledCounts, outputFile, 'Sheet', 'PooledCounts');
        writetable(ByContextCounts, outputFile, 'Sheet', 'ByContextCounts');
        writetable(IncludedValues, outputFile, 'Sheet', 'IncludedValues');
        fprintf('Saved main-room context place-field count workbook: %s\n', outputFile);
    catch ME
        warning('main-room context place-field count Excel export failed: %s', ME.message);
    end
end

function SummaryTable = compare_main_control_four_condition_aggregated_feature_sets(T, className, classOutDir, aggregationMode)

    [swrDefs, generalDefs] = build_feature_metric_defs(className);
    swrDefs = filter_metric_defs_by_available_data(T, swrDefs);
    generalDefs = filter_metric_defs_by_available_data(T, generalDefs);
    SummaryTable = table;
    roomOutDir = fullfile(classOutDir, 'room_comparisons');
    [aggregationLabel, fileTag, subsetLabel] = aggregation_mode_labels(aggregationMode);

    if ~exist(roomOutDir, 'dir')
        mkdir(roomOutDir);
    end

    figSWR = create_metric_tiled_figure(size(swrDefs, 1));
    for i = 1:size(swrDefs, 1)
        nexttile
        statsRows = plot_four_condition_aggregated_metric_panel(gca, T, swrDefs{i,1}, swrDefs{i,2}, swrDefs{i,3}, aggregationMode, i == 1);
        SummaryTable = [SummaryTable; statsRows]; %#ok<AGROW>
    end
    sgtitle(sprintf('%s cells: main-room/control cfos vs non-fos (%s)', className, aggregationLabel), 'Interpreter', 'none');
    save_figure_png(figSWR, fullfile(roomOutDir, sprintf('%s_SALT_MainControl_FourCondition_%s_SWRFeatures.png', className, fileTag)));
    close(figSWR);

    figGeneral = create_metric_tiled_figure(size(generalDefs, 1));
    for i = 1:size(generalDefs, 1)
        nexttile
        statsRows = plot_four_condition_aggregated_metric_panel(gca, T, generalDefs{i,1}, generalDefs{i,2}, generalDefs{i,3}, aggregationMode, i == 1);
        SummaryTable = [SummaryTable; statsRows]; %#ok<AGROW>
    end
    sgtitle(sprintf('%s cells: main-room/control cfos vs non-fos (%s)', className, aggregationLabel), 'Interpreter', 'none');
    save_figure_png(figGeneral, fullfile(roomOutDir, sprintf('%s_SALT_MainControl_FourCondition_%s_OtherFeatures.png', className, fileTag)));
    close(figGeneral);

    if ~isempty(SummaryTable)
        SummaryTable.CellClass = repmat(string(className), height(SummaryTable), 1);
        SummaryTable.ComparisonType = repmat(string(sprintf('main_control_four_condition_%s', aggregationMode)), height(SummaryTable), 1);
        SummaryTable.Subset = repmat(string(subsetLabel), height(SummaryTable), 1);
        SummaryTable.ControlRoom = repmat("All Control Rooms", height(SummaryTable), 1);
        SummaryTable = movevars(SummaryTable, {'CellClass','ComparisonType','Subset','ControlRoom'}, 'Before', 1);
        writetable(SummaryTable, fullfile(roomOutDir, sprintf('%s_SALT_MainControl_FourCondition_%s_FeatureSummary.csv', className, fileTag)));
    end
end

function mask = resolve_class_export_mask(T, className)

    if strcmpi(className, 'place cells')
        mask = T.CellClass == "principal" & T.IsPlaceCell;
    elseif strcmpi(className, 'non-place principal')
        mask = T.CellClass == "principal" & ~T.IsPlaceCell;
    else
        mask = T.CellClass == string(className);
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

function folderName = stage_b_export_folder_name(stageB_settings)

    stageB_settings = ensure_stage_b_settings_defaults(stageB_settings);
    folderName = sprintf('fixedLatency_%sms', sanitize_threshold_text(stageB_settings.fixed_latency_threshold_ms));
end

function T = ensure_numeric_table_columns(T, columnNames)

    for iColumn = 1:numel(columnNames)
        columnName = columnNames{iColumn};
        if ~ismember(columnName, T.Properties.VariableNames)
            T.(columnName) = nan(height(T), 1);
        end
    end
end

function values = pad_numeric_vector(rawValues, nRows)

    values = nan(nRows, 1);
    if isempty(rawValues)
        return
    end

    rawValues = double(rawValues(:));
    nCopy = min(numel(rawValues), nRows);
    values(1:nCopy) = rawValues(1:nCopy);
end

function T = assign_metric_by_cell_index(T, columnName, sessionRows, cellIdx, values)

    sessionRows = sessionRows(:);
    cellIdx = cellIdx(:);
    validIdx = isfinite(cellIdx) & cellIdx >= 1 & cellIdx <= numel(values);
    T.(columnName)(sessionRows(validIdx)) = values(cellIdx(validIdx));
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

function [swrDefs, generalDefs] = build_feature_metric_defs(className)

    if nargin < 1
        className = '';
    end

    swrDefs = { ...
        'S1_SSMI', 'S1', 'SSMI'; ...
        'S1_SFI',  'S1', 'SFI'; ...
        'S1_SPP',  'S1', 'SPP'; ...
        'S1_SMRpR', 'S1', 'Spikes per participated ripple'; ...
        'S2_SSMI', 'S2', 'SSMI'; ...
        'S2_SFI',  'S2', 'SFI'; ...
        'S2_SPP',  'S2', 'SPP'; ...
        'S2_SMRpR', 'S2', 'Spikes per participated ripple'};

    generalDefs = { ...
        'BurstIndex', 'Overall', 'Burst index'; ...
        'SpeedScore', 'Overall', 'Speed score'; ...
        'TemporalPeakRate', 'Overall', 'Temporal peak firing rate'; ...
        'AverageRate', 'Overall', 'Average firing rate'; ...
        'PlaceFieldSize', 'Overall', 'Place field size'; ...
        'SpatialInfo_A',     'Overall', 'Spatial information (bits/spike) Context A'; ... 
        'SpatialPeakRate', 'Overall', 'Spatial peak firing rate'; ...
        'AverageRate_OF1', 'OF1', 'Average firing rate'; ...
        'AverageRate_OF2', 'OF2', 'Average firing rate'; ...
        'AverageRate_OF3', 'OF3', 'Average firing rate'; ...
        'TemporalPeakRate_OF1', 'OF1', 'Temporal peak firing rate'; ...
        'TemporalPeakRate_OF2', 'OF2', 'Temporal peak firing rate'; ...
        'TemporalPeakRate_OF3', 'OF3', 'Temporal peak firing rate'; ...
        'SpatialPeakRate_OF1', 'OF1', 'Spatial peak firing rate'; ...
        'SpatialPeakRate_OF2', 'OF2', 'Spatial peak firing rate'; ...
        'SpatialPeakRate_OF3', 'OF3', 'Spatial peak firing rate'; ...
        'ThetaModulationScore', 'Overall', 'Theta modulation score'; ...
        'ThetaModulationTanaka', 'Overall', 'Theta modulation score (Tanaka)'; ...
        'IntrinsicFrequency', 'Overall', 'Intrinsic frequency'; ...
        'ThetaPhaseR', 'Overall', 'Theta phase-locking R'; ...
        'ThetaPhaseP', 'Overall', 'Theta phase-locking Rayleigh p'; ...
        'ThetaPhaseZ', 'Overall', 'Theta phase-locking Rayleigh z'; ...
        'ThetaPhaseNSpikes', 'Overall', 'Theta phase spike count'};

    if strcmpi(className, 'non-place principal')
        generalDefs = { ...
            'AverageRate_OFMean', 'Overall', 'Average firing rate (mean OF1-3)'; ...
            'TemporalPeakRate', 'Overall', 'Temporal peak firing rate'; ...
            'SpatialPeakRate_OFMax', 'Overall', 'Spatial peak firing rate (max OF1-3)'; ...
            'AverageRate_OF1', 'OF1', 'Average firing rate'; ...
            'AverageRate_OF2', 'OF2', 'Average firing rate'; ...
            'AverageRate_OF3', 'OF3', 'Average firing rate'; ...
            'TemporalPeakRate_OF1', 'OF1', 'Temporal peak firing rate'; ...
            'TemporalPeakRate_OF2', 'OF2', 'Temporal peak firing rate'; ...
            'TemporalPeakRate_OF3', 'OF3', 'Temporal peak firing rate'; ...
            'SpatialPeakRate_OF1', 'OF1', 'Spatial peak firing rate'; ...
            'SpatialPeakRate_OF2', 'OF2', 'Spatial peak firing rate'; ...
            'SpatialPeakRate_OF3', 'OF3', 'Spatial peak firing rate'; ...
            'BurstIndex_OFMean', 'Overall', 'Burst index (mean OF1-3)'; ...
            'ThetaModulationScore_OFMean', 'Overall', 'Theta modulation score (mean OF1-3)'; ...
            'ThetaModulationTanaka_OFMean', 'Overall', 'Theta modulation score (Tanaka mean OF1-3)'; ...
            'IntrinsicFrequency_OFMean', 'Overall', 'Intrinsic frequency (mean OF1-3)'; ...
            'ThetaPhaseR', 'Overall', 'Theta phase-locking R'; ...
            'ThetaPhaseP', 'Overall', 'Theta phase-locking Rayleigh p'; ...
            'ThetaPhaseZ', 'Overall', 'Theta phase-locking Rayleigh z'; ...
            'ThetaPhaseNSpikes', 'Overall', 'Theta phase spike count'; ...
            'SpeedScore', 'Overall', 'Speed score'; ...
            'SpatialInfo_A', 'Overall', 'Spatial information (bits/spike) Context A'};
    elseif strcmpi(className, 'place cells')
        generalDefs = [generalDefs; { ...
            'PlaceFieldNumber', 'Overall', 'Place field number'; ...
            'PlaceFieldSize_OF1', 'OF1', 'Place field size'; ...
            'PlaceFieldSize_OF2', 'OF2', 'Place field size'; ...
            'PlaceFieldSize_OF3', 'OF3', 'Place field size'; ...
            'PlaceFieldNumber_OF1', 'OF1', 'Place field number'; ...
            'PlaceFieldNumber_OF2', 'OF2', 'Place field number'; ...
            'PlaceFieldNumber_OF3', 'OF3', 'Place field number'; ...
            'ZCoherence_OF1', 'OF1', 'Place field z-coherence'; ...
            'ZCoherence_OF2', 'OF2', 'Place field z-coherence'; ...
            'ZCoherence_OF3', 'OF3', 'Place field z-coherence'; ...
            'Sparseness_OF1', 'OF1', 'Place field sparseness'; ...
            'Sparseness_OF2', 'OF2', 'Place field sparseness'; ...
            'Sparseness_OF3', 'OF3', 'Place field sparseness'; ...
            'Selectivity_OF1', 'OF1', 'Place field selectivity'; ...
            'Selectivity_OF2', 'OF2', 'Place field selectivity'; ...
            'Selectivity_OF3', 'OF3', 'Place field selectivity'; ...
            'SpatialInfo_OF1', 'OF1', 'Spatial information (bits/spike)'; ...
            'SpatialInfo_OF2', 'OF2', 'Spatial information (bits/spike)'; ...
            'SpatialInfo_OF3', 'OF3', 'Spatial information (bits/spike)'; ...
            'SpatialInfoRate_OF1', 'OF1', 'Spatial information rate (bits/sec)'; ...
            'SpatialInfoRate_OF2', 'OF2', 'Spatial information rate (bits/sec)'; ...
            'SpatialInfoRate_OF3', 'OF3', 'Spatial information rate (bits/sec)'; ...
            'Remapping', 'Overall', 'Remapping (OF1 vs OF2)'; ...
            'Stability', 'Overall', 'Stability (OF1 vs OF3)'}];
    end
end

function defs = filter_metric_defs_by_available_data(T, defs)

    if isempty(defs)
        return
    end

    keepMask = false(size(defs, 1), 1);
    variableNames = string(T.Properties.VariableNames);

    for iDef = 1:size(defs, 1)
        fieldName = string(defs{iDef,1});
        if ~any(variableNames == fieldName)
            continue
        end

        values = T.(fieldName);
        keepMask(iDef) = any(isfinite(values));
    end

    defs = defs(keepMask, :);
end

function statsRow = plot_direct_vs_rest_metric_panel(ax, T, fieldName, stateName, featureName, showLegend)

    values = T.(fieldName);
    classifiedMask = T.FinalLabel ~= "unclassified";
    directMask = classifiedMask & T.FinalLabel == "direct" & ~isnan(values);
    indirectMask = classifiedMask & T.FinalLabel == "indirect" & ~isnan(values);
    untaggedMask = classifiedMask & T.FinalLabel == "untagged" & ~isnan(values);
    restMask = indirectMask | untaggedMask;

    xRestUntagged = 1 + 0.16 * (rand(nnz(untaggedMask), 1) - 0.5);
    xRestIndirect = 1 + 0.16 * (rand(nnz(indirectMask), 1) - 0.5);
    xDirect = 2 + 0.16 * (rand(nnz(directMask), 1) - 0.5);

    hold(ax, 'on');

    hUntagged = gobjects(0);
    hIndirect = gobjects(0);
    hDirect = gobjects(0);

    if any(untaggedMask)
        hUntagged = scatter(ax, xRestUntagged, values(untaggedMask), 28, [0.55 0.55 0.55], 'filled', ...
            'MarkerFaceAlpha', 0.7, 'MarkerEdgeColor', 'none');
    end

    if any(indirectMask)
        hIndirect = scatter(ax, xRestIndirect, values(indirectMask), 34, [0.80 0.80 0.80], 'filled', ...
            'MarkerFaceAlpha', 0.8, 'MarkerEdgeColor', [0.15 0.55 0.15], 'LineWidth', 1.0);
    end

    if any(restMask)
        plot_group_summary(ax, 1, values(restMask), [0.25 0.25 0.25]);
    end

    if any(directMask)
        hDirect = scatter(ax, xDirect, values(directMask), 34, [0.2 0.7 0.2], 'filled', ...
            'MarkerFaceAlpha', 0.8, 'MarkerEdgeColor', [0.1 0.1 0.1], 'LineWidth', 0.3);
        plot_group_summary(ax, 2, values(directMask), [0.1 0.45 0.1]);
    end

    [pVal, zVal] = run_ranksum_safe(values(restMask), values(directMask));

    set(ax, 'XLim', [0.5 2.5], 'XTick', [1 2], 'XTickLabel', {'untagged + indirect','direct'});
    ylabel(ax, featureName);
    title(ax, sprintf('%s %s\np = %s', stateName, featureName, format_p_value(pVal)), 'Interpreter', 'none');
    grid(ax, 'on');
    box(ax, 'on');

    if showLegend
        legendHandles = gobjects(0);
        legendLabels = {};
        if ~isempty(hUntagged)
            legendHandles(end + 1) = hUntagged; %#ok<AGROW>
            legendLabels{end + 1} = 'untagged'; %#ok<AGROW>
        end
        if ~isempty(hIndirect)
            legendHandles(end + 1) = hIndirect; %#ok<AGROW>
            legendLabels{end + 1} = 'indirect'; %#ok<AGROW>
        end
        if ~isempty(hDirect)
            legendHandles(end + 1) = hDirect; %#ok<AGROW>
            legendLabels{end + 1} = 'direct'; %#ok<AGROW>
        end
        if ~isempty(legendHandles)
            legend(ax, legendHandles, legendLabels, 'Location', 'best');
        end
    end

    statsRow = build_comparison_row(stateName, featureName, ...
        "untagged + indirect", "direct", values(restMask), values(directMask), nnz(indirectMask), pVal, zVal);
end

function figHandle = create_metric_tiled_figure(nPanels)

    nPanels = max(1, double(nPanels));
    nCols = min(3, nPanels);
    nRows = ceil(nPanels / nCols);

    figWidth = 420 * nCols + 80;
    figHeight = 300 * nRows + 120;

    figHandle = figure('Color', 'w', 'Position', [100 100 figWidth figHeight]);
    tiledlayout(nRows, nCols, 'Padding', 'compact', 'TileSpacing', 'compact');
end

function save_figure_png(figHandle, filename)
    [folder,~,~] = fileparts(filename);
    if ~exist(folder, 'dir')
        mkdir(folder);
    end
    exportgraphics(figHandle, filename, 'Resolution', 300);
end

function statsRow = plot_room_group_metric_panel(ax, T, fieldName, stateName, featureName, labelName)

    values = T.(fieldName);
    labelMask = resolve_label_mask(T, labelName) & ~isnan(values);
    mainMask = labelMask & T.RoomGroup == "main";
    controlMask = labelMask & T.RoomGroup == "control";
    highlightIndirect = strcmp(labelName, 'non_direct') || strcmp(labelName, 'rest');
    indirectMainMask = highlightIndirect & mainMask & T.FinalLabel == "indirect";
    indirectControlMask = highlightIndirect & controlMask & T.FinalLabel == "indirect";
    untaggedMainMask = mainMask & ~indirectMainMask;
    untaggedControlMask = controlMask & ~indirectControlMask;

    hold(ax, 'on');

    if any(untaggedMainMask)
        xMainUntagged = 1 + 0.16 * (rand(nnz(untaggedMainMask), 1) - 0.5);
        scatter(ax, xMainUntagged, values(untaggedMainMask), 30, [0.55 0.55 0.55], 'filled', ...
            'MarkerFaceAlpha', 0.75, 'MarkerEdgeColor', [0.1 0.1 0.1], 'LineWidth', 0.8);
    end
    if any(indirectMainMask)
        xMainIndirect = 1 + 0.16 * (rand(nnz(indirectMainMask), 1) - 0.5);
        scatter(ax, xMainIndirect, values(indirectMainMask), 30, [0.55 0.55 0.55], 'filled', ...
            'MarkerFaceAlpha', 0.75, 'MarkerEdgeColor', [0.15 0.55 0.15], 'LineWidth', 1.0);
    end
    if any(mainMask)
        plot_group_summary(ax, 1, values(mainMask), [0.25 0.25 0.25]);
    end

    if any(untaggedControlMask)
        xControlUntagged = 2 + 0.16 * (rand(nnz(untaggedControlMask), 1) - 0.5);
        scatter(ax, xControlUntagged, values(untaggedControlMask), 30, [0.90 0.60 0.20], 'filled', ...
            'MarkerFaceAlpha', 0.80, 'MarkerEdgeColor', [0.1 0.1 0.1], 'LineWidth', 0.8);
    end
    if any(indirectControlMask)
        xControlIndirect = 2 + 0.16 * (rand(nnz(indirectControlMask), 1) - 0.5);
        scatter(ax, xControlIndirect, values(indirectControlMask), 30, [0.90 0.60 0.20], 'filled', ...
            'MarkerFaceAlpha', 0.80, 'MarkerEdgeColor', [0.15 0.55 0.15], 'LineWidth', 1.0);
    end
    if any(controlMask)
        plot_group_summary(ax, 2, values(controlMask), [0.75 0.35 0.10]);
    end

    [pVal, zVal] = run_ranksum_safe(values(mainMask), values(controlMask));

    set(ax, 'XLim', [0.5 2.5], 'XTick', [1 2], 'XTickLabel', {'main-room','control rooms'});
    ylabel(ax, featureName);
    title(ax, sprintf('%s %s\np = %s', stateName, featureName, format_p_value(pVal)), 'Interpreter', 'none');
    grid(ax, 'on');
    box(ax, 'on');

    statsRow = build_comparison_row(stateName, featureName, ...
        "main-room", "control rooms", values(mainMask), values(controlMask), ...
        nnz(indirectMainMask) + nnz(indirectControlMask), pVal, zVal);
end

function labelMask = resolve_label_mask(T, labelName)

    if strcmp(labelName, 'tagged')
        labelMask = T.StageA_Label == "tagged";
    elseif strcmp(labelName, 'non_direct') || strcmp(labelName, 'rest')
        labelMask = T.FinalLabel == "untagged" | T.FinalLabel == "indirect";
    else
        labelMask = T.FinalLabel == labelName;
    end
end

function labelText = label_display_name(labelName)

    if strcmp(labelName, 'non_direct') || strcmp(labelName, 'rest')
        labelText = 'untagged + indirect';
    else
        labelText = labelName;
    end
end

function statsRow = plot_specific_room_metric_panel(ax, T, fieldName, stateName, featureName, labelName, controlRoomId, controlRoomDisplayName)

    values = T.(fieldName);
    labelMask = resolve_label_mask(T, labelName) & ~isnan(values);
    mainMask = labelMask & T.RoomID == 15;
    controlMask = labelMask & T.RoomID == controlRoomId;
    highlightIndirect = strcmp(labelName, 'non_direct') || strcmp(labelName, 'rest');
    indirectMainMask = highlightIndirect & mainMask & T.FinalLabel == "indirect";
    indirectControlMask = highlightIndirect & controlMask & T.FinalLabel == "indirect";
    untaggedMainMask = mainMask & ~indirectMainMask;
    untaggedControlMask = controlMask & ~indirectControlMask;

    hold(ax, 'on');

    if any(untaggedMainMask)
        xMainUntagged = 1 + 0.16 * (rand(nnz(untaggedMainMask), 1) - 0.5);
        scatter(ax, xMainUntagged, values(untaggedMainMask), 30, [0.55 0.55 0.55], 'filled', ...
            'MarkerFaceAlpha', 0.75, 'MarkerEdgeColor', [0.1 0.1 0.1], 'LineWidth', 0.8);
    end
    if any(indirectMainMask)
        xMainIndirect = 1 + 0.16 * (rand(nnz(indirectMainMask), 1) - 0.5);
        scatter(ax, xMainIndirect, values(indirectMainMask), 30, [0.55 0.55 0.55], 'filled', ...
            'MarkerFaceAlpha', 0.75, 'MarkerEdgeColor', [0.15 0.55 0.15], 'LineWidth', 1.0);
    end
    if any(mainMask)
        plot_group_summary(ax, 1, values(mainMask), [0.25 0.25 0.25]);
    end

    if any(untaggedControlMask)
        xControlUntagged = 2 + 0.16 * (rand(nnz(untaggedControlMask), 1) - 0.5);
        scatter(ax, xControlUntagged, values(untaggedControlMask), 30, [0.20 0.45 0.85], 'filled', ...
            'MarkerFaceAlpha', 0.80, 'MarkerEdgeColor', [0.1 0.1 0.1], 'LineWidth', 0.8);
    end
    if any(indirectControlMask)
        xControlIndirect = 2 + 0.16 * (rand(nnz(indirectControlMask), 1) - 0.5);
        scatter(ax, xControlIndirect, values(indirectControlMask), 30, [0.20 0.45 0.85], 'filled', ...
            'MarkerFaceAlpha', 0.80, 'MarkerEdgeColor', [0.15 0.55 0.15], 'LineWidth', 1.0);
    end
    if any(controlMask)
        plot_group_summary(ax, 2, values(controlMask), [0.10 0.30 0.65]);
    end

    [pVal, zVal] = run_ranksum_safe(values(mainMask), values(controlMask));

    set(ax, 'XLim', [0.5 2.5], 'XTick', [1 2], 'XTickLabel', {'main-room', controlRoomDisplayName});
    ylabel(ax, featureName);
    title(ax, sprintf('%s %s\np = %s', stateName, featureName, format_p_value(pVal)), 'Interpreter', 'none');
    grid(ax, 'on');
    box(ax, 'on');

    statsRow = build_comparison_row(stateName, featureName, ...
        "main-room", string(controlRoomDisplayName), values(mainMask), values(controlMask), ...
        nnz(indirectMainMask) + nnz(indirectControlMask), pVal, zVal);
end

function export_four_condition_csv_columns(T, fieldName, stateName, featureName, className, csvExportDir)

    values = T.(fieldName);
    classifiedMask = T.FinalLabel ~= "unclassified";

    mainDirectMask = classifiedMask & T.RoomGroup == "main" & T.FinalLabel == "direct";
    mainRestMask = classifiedMask & T.RoomGroup == "main" & ...
        (T.FinalLabel == "untagged" | T.FinalLabel == "indirect");
    controlDirectMask = classifiedMask & T.RoomGroup == "control" & T.FinalLabel == "direct";
    controlRestMask = classifiedMask & T.RoomGroup == "control" & ...
        (T.FinalLabel == "untagged" | T.FinalLabel == "indirect");

    classFolder = four_condition_class_folder_name(className);
    graphFolder = sanitize_filename(sprintf('%s_%s', stateName, featureName));
    exportDir = fullfile(csvExportDir, classFolder, graphFolder);

    if ~exist(exportDir, 'dir')
        mkdir(exportDir);
    end

    write_column_csv(values(mainRestMask), fullfile(exportDir, 'main_room_non_fos.csv'));
    write_column_csv(values(mainDirectMask), fullfile(exportDir, 'main_room_cfos.csv'));
    write_column_csv(values(controlRestMask), fullfile(exportDir, 'control_non_fos.csv'));
    write_column_csv(values(controlDirectMask), fullfile(exportDir, 'control_cfos.csv'));
end

function statsRows = plot_four_condition_metric_panel(ax, T, fieldName, stateName, featureName, showLegend)

    values = T.(fieldName);
    classifiedMask = T.FinalLabel ~= "unclassified" & ~isnan(values);

    mainDirectMask = classifiedMask & T.RoomGroup == "main" & T.FinalLabel == "direct";
    mainIndirectMask = classifiedMask & T.RoomGroup == "main" & T.FinalLabel == "indirect";
    mainUntaggedMask = classifiedMask & T.RoomGroup == "main" & T.FinalLabel == "untagged";
    controlDirectMask = classifiedMask & T.RoomGroup == "control" & T.FinalLabel == "direct";
    controlIndirectMask = classifiedMask & T.RoomGroup == "control" & T.FinalLabel == "indirect";
    controlUntaggedMask = classifiedMask & T.RoomGroup == "control" & T.FinalLabel == "untagged";

    mainRestMask = mainIndirectMask | mainUntaggedMask;
    controlRestMask = controlIndirectMask | controlUntaggedMask;

    hold(ax, 'on');

    hMainRest = gobjects(0);
    hMainDirect = gobjects(0);
    hControlRest = gobjects(0);
    hControlDirect = gobjects(0);

    if any(mainUntaggedMask)
        x = 1 + 0.16 * (rand(nnz(mainUntaggedMask), 1) - 0.5);
        hMainRest = scatter(ax, x, values(mainUntaggedMask), 28, [0.55 0.55 0.55], 'filled', ...
            'MarkerFaceAlpha', 0.75, 'MarkerEdgeColor', [0.1 0.1 0.1], 'LineWidth', 0.8);
    end
    if any(mainIndirectMask)
        x = 1 + 0.16 * (rand(nnz(mainIndirectMask), 1) - 0.5);
        scatter(ax, x, values(mainIndirectMask), 30, [0.55 0.55 0.55], 'filled', ...
            'MarkerFaceAlpha', 0.75, 'MarkerEdgeColor', [0.15 0.55 0.15], 'LineWidth', 1.0);
    end
    if any(mainRestMask)
        plot_group_summary(ax, 1, values(mainRestMask), [0.25 0.25 0.25]);
    end

    if any(mainDirectMask)
        x = 2 + 0.16 * (rand(nnz(mainDirectMask), 1) - 0.5);
        hMainDirect = scatter(ax, x, values(mainDirectMask), 34, [0.2 0.7 0.2], 'filled', ...
            'MarkerFaceAlpha', 0.8, 'MarkerEdgeColor', [0.1 0.1 0.1], 'LineWidth', 0.3);
        plot_group_summary(ax, 2, values(mainDirectMask), [0.1 0.45 0.1]);
    end

    if any(controlUntaggedMask)
        x = 3 + 0.16 * (rand(nnz(controlUntaggedMask), 1) - 0.5);
        hControlRest = scatter(ax, x, values(controlUntaggedMask), 28, [0.90 0.60 0.20], 'filled', ...
            'MarkerFaceAlpha', 0.80, 'MarkerEdgeColor', [0.1 0.1 0.1], 'LineWidth', 0.8);
    end
    if any(controlIndirectMask)
        x = 3 + 0.16 * (rand(nnz(controlIndirectMask), 1) - 0.5);
        scatter(ax, x, values(controlIndirectMask), 30, [0.90 0.60 0.20], 'filled', ...
            'MarkerFaceAlpha', 0.80, 'MarkerEdgeColor', [0.15 0.55 0.15], 'LineWidth', 1.0);
    end
    if any(controlRestMask)
        plot_group_summary(ax, 3, values(controlRestMask), [0.75 0.35 0.10]);
    end

    if any(controlDirectMask)
        x = 4 + 0.16 * (rand(nnz(controlDirectMask), 1) - 0.5);
        hControlDirect = scatter(ax, x, values(controlDirectMask), 34, [0.20 0.45 0.85], 'filled', ...
            'MarkerFaceAlpha', 0.85, 'MarkerEdgeColor', [0.1 0.1 0.1], 'LineWidth', 0.3);
        plot_group_summary(ax, 4, values(controlDirectMask), [0.10 0.30 0.65]);
    end

    set(ax, 'XLim', [0.5 4.5], 'XTick', 1:4, ...
        'XTickLabel', {'main-room non-fos','main-room cfos','control non-fos','control cfos'});
    xtickangle(ax, 20);
    ylabel(ax, featureName);
    title(ax, sprintf('%s %s', stateName, featureName), 'Interpreter', 'none');
    grid(ax, 'on');
    box(ax, 'on');

    [pMain, zMain] = run_ranksum_safe(values(mainRestMask), values(mainDirectMask));
    [pControl, zControl] = run_ranksum_safe(values(controlRestMask), values(controlDirectMask));
    [pDirectRoom, zDirectRoom] = run_ranksum_safe(values(mainDirectMask), values(controlDirectMask));
    [pRestRoom, zRestRoom] = run_ranksum_safe(values(mainRestMask), values(controlRestMask));

    yLimits = ylim(ax);
    xText = 0.02;
    yText = 0.98;
    text(ax, xText, yText, sprintf(['main-room non-fos vs cfos: %s\n' ...
        'control non-fos vs cfos: %s\n' ...
        'main-room cfos vs control cfos: %s\n' ...
        'main-room non-fos vs control non-fos: %s'], ...
        format_p_value(pMain), format_p_value(pControl), ...
        format_p_value(pDirectRoom), format_p_value(pRestRoom)), ...
        'Units', 'normalized', 'HorizontalAlignment', 'left', 'VerticalAlignment', 'top', ...
        'FontSize', 8, 'BackgroundColor', 'w', 'Margin', 4);

    if showLegend
        legendHandles = gobjects(0);
        legendLabels = {};
        if ~isempty(hMainRest)
            legendHandles(end + 1) = hMainRest; %#ok<AGROW>
            legendLabels{end + 1} = 'main-room non-fos'; %#ok<AGROW>
        end
        if any(mainIndirectMask)
            hIndirect = scatter(ax, nan, nan, 30, [0.75 0.75 0.75], 'filled', ...
                'MarkerEdgeColor', [0.15 0.55 0.15], 'LineWidth', 1.0);
            legendHandles(end + 1) = hIndirect; %#ok<AGROW>
            legendLabels{end + 1} = 'indirect'; %#ok<AGROW>
        end
        if ~isempty(hMainDirect)
            legendHandles(end + 1) = hMainDirect; %#ok<AGROW>
            legendLabels{end + 1} = 'main-room cfos'; %#ok<AGROW>
        end
        if ~isempty(hControlRest)
            legendHandles(end + 1) = hControlRest; %#ok<AGROW>
            legendLabels{end + 1} = 'control non-fos'; %#ok<AGROW>
        end
        if ~isempty(hControlDirect)
            legendHandles(end + 1) = hControlDirect; %#ok<AGROW>
            legendLabels{end + 1} = 'control cfos'; %#ok<AGROW>
        end
        if ~isempty(legendHandles)
            legend(ax, legendHandles, legendLabels, 'Location', 'best');
        end
    end

    statsRows = table;
    statsRows = [statsRows; build_comparison_row(stateName, featureName, ...
        "main-room non-fos", "main-room cfos", values(mainRestMask), values(mainDirectMask), ...
        nnz(mainIndirectMask), pMain, zMain)]; %#ok<AGROW>
    statsRows = [statsRows; build_comparison_row(stateName, featureName, ...
        "control non-fos", "control cfos", values(controlRestMask), values(controlDirectMask), ...
        nnz(controlIndirectMask), pControl, zControl)]; %#ok<AGROW>
    statsRows = [statsRows; build_comparison_row(stateName, featureName, ...
        "main-room cfos", "control cfos", values(mainDirectMask), values(controlDirectMask), ...
        0, pDirectRoom, zDirectRoom)]; %#ok<AGROW>
    statsRows = [statsRows; build_comparison_row(stateName, featureName, ...
        "main-room non-fos", "control non-fos", values(mainRestMask), values(controlRestMask), ...
        nnz(mainIndirectMask) + nnz(controlIndirectMask), pRestRoom, zRestRoom)]; %#ok<AGROW>
end

function [IncludedValues, contextNames, contextLabels] = build_main_room_context_place_field_values(T)

    contextNames = ["OF1"; "OF2"; "OF3"];
    contextLabels = ["Context A"; "Context B"; "Context A revisit"];
    pfColumns = {'PlaceFieldNumber_OF1', 'PlaceFieldNumber_OF2', 'PlaceFieldNumber_OF3'};
    if ismember('RoomName', T.Properties.VariableNames)
        roomName = T.RoomName;
    else
        roomName = repmat("main-room", height(T), 1);
    end

    IncludedValues = table( ...
        nan(0,1), nan(0,1), nan(0,1), strings(0,1), strings(0,1), ...
        strings(0,1), strings(0,1), strings(0,1), strings(0,1), strings(0,1), ...
        nan(0,1), strings(0,1), ...
        'VariableNames', {'SessionIndex','INumber','Cell','AnimalID','RoomName', ...
        'CellClass','FinalLabel','ExportGroup','Context','ContextLabel', ...
        'PlaceFieldNumber','PlaceFieldBin'});

    exportGroup = main_room_context_place_field_export_group(T.FinalLabel);
    baseMask = T.RoomID == 15 & T.CellClass == "principal" & exportGroup ~= "unclassified";

    for iContext = 1:numel(contextNames)
        fieldCounts = T.(pfColumns{iContext});
        contextMask = baseMask & isfinite(fieldCounts) & fieldCounts > 0;
        if ~any(contextMask)
            continue
        end

        nRows = nnz(contextMask);
        contextCounts = fieldCounts(contextMask);
        contextTable = table( ...
            T.SessionIndex(contextMask), ...
            T.INumber(contextMask), ...
            T.Cell(contextMask), ...
            T.AnimalID(contextMask), ...
            roomName(contextMask), ...
            T.CellClass(contextMask), ...
            T.FinalLabel(contextMask), ...
            exportGroup(contextMask), ...
            repmat(contextNames(iContext), nRows, 1), ...
            repmat(contextLabels(iContext), nRows, 1), ...
            contextCounts(:), ...
            main_room_context_place_field_bin(contextCounts), ...
            'VariableNames', IncludedValues.Properties.VariableNames);
        IncludedValues = [IncludedValues; contextTable]; %#ok<AGROW>
    end
end

function Summary = summarize_main_room_context_place_field_values(IncludedValues, contextNames, contextLabels, splitByContext)

    exportGroups = ["main_room_non_fos"; "main_room_cfos"];
    saltGroups = ["untagged + indirect"; "direct"];

    RoomName = strings(0,1);
    ExportGroup = strings(0,1);
    SALTGroup = strings(0,1);
    Context = strings(0,1);
    ContextLabel = strings(0,1);
    NContextValues = [];
    NUniqueCells = [];
    Count_1PlaceField = [];
    Count_2PlaceFields = [];
    Count_3PlusPlaceFields = [];
    Fraction_1PlaceField = [];
    Fraction_2PlaceFields = [];
    Fraction_3PlusPlaceFields = [];

    for iGroup = 1:numel(exportGroups)
        if splitByContext
            for iContext = 1:numel(contextNames)
                valueMask = IncludedValues.ExportGroup == exportGroups(iGroup) & ...
                    IncludedValues.Context == contextNames(iContext);
                [nValues, nCells, nOne, nTwo, nThreePlus, fOne, fTwo, fThreePlus] = ...
                    count_main_room_context_place_field_bins(IncludedValues(valueMask, :));

                RoomName = [RoomName; "main-room"]; %#ok<AGROW>
                ExportGroup = [ExportGroup; exportGroups(iGroup)]; %#ok<AGROW>
                SALTGroup = [SALTGroup; saltGroups(iGroup)]; %#ok<AGROW>
                Context = [Context; contextNames(iContext)]; %#ok<AGROW>
                ContextLabel = [ContextLabel; contextLabels(iContext)]; %#ok<AGROW>
                NContextValues = [NContextValues; nValues]; %#ok<AGROW>
                NUniqueCells = [NUniqueCells; nCells]; %#ok<AGROW>
                Count_1PlaceField = [Count_1PlaceField; nOne]; %#ok<AGROW>
                Count_2PlaceFields = [Count_2PlaceFields; nTwo]; %#ok<AGROW>
                Count_3PlusPlaceFields = [Count_3PlusPlaceFields; nThreePlus]; %#ok<AGROW>
                Fraction_1PlaceField = [Fraction_1PlaceField; fOne]; %#ok<AGROW>
                Fraction_2PlaceFields = [Fraction_2PlaceFields; fTwo]; %#ok<AGROW>
                Fraction_3PlusPlaceFields = [Fraction_3PlusPlaceFields; fThreePlus]; %#ok<AGROW>
            end
        else
            valueMask = IncludedValues.ExportGroup == exportGroups(iGroup);
            [nValues, nCells, nOne, nTwo, nThreePlus, fOne, fTwo, fThreePlus] = ...
                count_main_room_context_place_field_bins(IncludedValues(valueMask, :));

            RoomName = [RoomName; "main-room"]; %#ok<AGROW>
            ExportGroup = [ExportGroup; exportGroups(iGroup)]; %#ok<AGROW>
            SALTGroup = [SALTGroup; saltGroups(iGroup)]; %#ok<AGROW>
            Context = [Context; "OF1-OF3 pooled"]; %#ok<AGROW>
            ContextLabel = [ContextLabel; "All contexts pooled"]; %#ok<AGROW>
            NContextValues = [NContextValues; nValues]; %#ok<AGROW>
            NUniqueCells = [NUniqueCells; nCells]; %#ok<AGROW>
            Count_1PlaceField = [Count_1PlaceField; nOne]; %#ok<AGROW>
            Count_2PlaceFields = [Count_2PlaceFields; nTwo]; %#ok<AGROW>
            Count_3PlusPlaceFields = [Count_3PlusPlaceFields; nThreePlus]; %#ok<AGROW>
            Fraction_1PlaceField = [Fraction_1PlaceField; fOne]; %#ok<AGROW>
            Fraction_2PlaceFields = [Fraction_2PlaceFields; fTwo]; %#ok<AGROW>
            Fraction_3PlusPlaceFields = [Fraction_3PlusPlaceFields; fThreePlus]; %#ok<AGROW>
        end
    end

    Summary = table(RoomName, ExportGroup, SALTGroup, Context, ContextLabel, ...
        NContextValues, NUniqueCells, Count_1PlaceField, Count_2PlaceFields, ...
        Count_3PlusPlaceFields, Fraction_1PlaceField, Fraction_2PlaceFields, ...
        Fraction_3PlusPlaceFields);
end

function Notes = main_room_context_place_field_notes_table()

    Item = [ ...
        "Selection"; ...
        "Pooled analysis"; ...
        "main_room_cfos"; ...
        "main_room_non_fos"; ...
        "Count columns"];
    Description = [ ...
        "main-room, CellClass == principal, classified SALT label, and context-specific place-field number > 0."; ...
        "OF1, OF2, and OF3 place-cell entries are stacked; one cell can contribute up to three context values."; ...
        "FinalLabel == direct."; ...
        "FinalLabel == untagged or indirect."; ...
        "Counts use the context-specific place-field number: 1, 2, or 3+ place fields."];
    Notes = table(Item, Description);
end

function statsRows = plot_four_condition_aggregated_metric_panel(ax, T, fieldName, stateName, featureName, aggregationMode, showLegend)

    values = T.(fieldName);
    classifiedMask = T.FinalLabel ~= "unclassified" & ~isnan(values);

    mainDirectMask = classifiedMask & T.RoomGroup == "main" & T.FinalLabel == "direct";
    mainIndirectMask = classifiedMask & T.RoomGroup == "main" & T.FinalLabel == "indirect";
    mainUntaggedMask = classifiedMask & T.RoomGroup == "main" & T.FinalLabel == "untagged";
    controlDirectMask = classifiedMask & T.RoomGroup == "control" & T.FinalLabel == "direct";
    controlIndirectMask = classifiedMask & T.RoomGroup == "control" & T.FinalLabel == "indirect";
    controlUntaggedMask = classifiedMask & T.RoomGroup == "control" & T.FinalLabel == "untagged";

    mainRestMask = mainIndirectMask | mainUntaggedMask;
    controlRestMask = controlIndirectMask | controlUntaggedMask;

    [mainRestValues, mainRestHasIndirect] = aggregate_group_means(T, values, mainRestMask, mainIndirectMask, aggregationMode);
    [mainDirectValues, ~] = aggregate_group_means(T, values, mainDirectMask, false(size(mainDirectMask)), aggregationMode);
    [controlRestValues, controlRestHasIndirect] = aggregate_group_means(T, values, controlRestMask, controlIndirectMask, aggregationMode);
    [controlDirectValues, ~] = aggregate_group_means(T, values, controlDirectMask, false(size(controlDirectMask)), aggregationMode);

    hold(ax, 'on');

    hMainRest = gobjects(0);
    hMainDirect = gobjects(0);
    hControlRest = gobjects(0);
    hControlDirect = gobjects(0);

    [hMainRest, hMainIndirect] = scatter_aggregated_group(ax, 1, mainRestValues, mainRestHasIndirect, [0.55 0.55 0.55]);
    if ~isempty(mainRestValues)
        plot_group_summary(ax, 1, mainRestValues, [0.25 0.25 0.25]);
    end

    if ~isempty(mainDirectValues)
        x = 2 + 0.16 * (rand(numel(mainDirectValues), 1) - 0.5);
        hMainDirect = scatter(ax, x, mainDirectValues, 36, [0.2 0.7 0.2], 'filled', ...
            'MarkerFaceAlpha', 0.85, 'MarkerEdgeColor', [0.1 0.1 0.1], 'LineWidth', 0.3);
        plot_group_summary(ax, 2, mainDirectValues, [0.1 0.45 0.1]);
    end

    [hControlRest, hControlIndirect] = scatter_aggregated_group(ax, 3, controlRestValues, controlRestHasIndirect, [0.90 0.60 0.20]);
    if ~isempty(controlRestValues)
        plot_group_summary(ax, 3, controlRestValues, [0.75 0.35 0.10]);
    end

    if ~isempty(controlDirectValues)
        x = 4 + 0.16 * (rand(numel(controlDirectValues), 1) - 0.5);
        hControlDirect = scatter(ax, x, controlDirectValues, 36, [0.20 0.45 0.85], 'filled', ...
            'MarkerFaceAlpha', 0.85, 'MarkerEdgeColor', [0.1 0.1 0.1], 'LineWidth', 0.3);
        plot_group_summary(ax, 4, controlDirectValues, [0.10 0.30 0.65]);
    end

    set(ax, 'XLim', [0.5 4.5], 'XTick', 1:4, ...
        'XTickLabel', {'main-room non-fos','main-room cfos','control non-fos','control cfos'});
    xtickangle(ax, 20);
    ylabel(ax, sprintf('%s mean', featureName));
    title(ax, sprintf('%s %s (%s)', stateName, featureName, aggregationMode), 'Interpreter', 'none');
    grid(ax, 'on');
    box(ax, 'on');

    [pMain, zMain] = run_ranksum_safe(mainRestValues, mainDirectValues);
    [pControl, zControl] = run_ranksum_safe(controlRestValues, controlDirectValues);
    [pDirectRoom, zDirectRoom] = run_ranksum_safe(mainDirectValues, controlDirectValues);
    [pRestRoom, zRestRoom] = run_ranksum_safe(mainRestValues, controlRestValues);

    text(ax, 0.02, 0.98, sprintf(['main-room non-fos vs cfos: %s\n' ...
        'control non-fos vs cfos: %s\n' ...
        'main-room cfos vs control cfos: %s\n' ...
        'main-room non-fos vs control non-fos: %s'], ...
        format_p_value(pMain), format_p_value(pControl), ...
        format_p_value(pDirectRoom), format_p_value(pRestRoom)), ...
        'Units', 'normalized', 'HorizontalAlignment', 'left', 'VerticalAlignment', 'top', ...
        'FontSize', 8, 'BackgroundColor', 'w', 'Margin', 4);

    if showLegend
        legendHandles = gobjects(0);
        legendLabels = {};
        if ~isempty(hMainRest)
            legendHandles(end + 1) = hMainRest; %#ok<AGROW>
            legendLabels{end + 1} = sprintf('main-room non-fos (%s mean)', aggregationMode); %#ok<AGROW>
        end
        if ~isempty(hMainIndirect)
            legendHandles(end + 1) = hMainIndirect; %#ok<AGROW>
            legendLabels{end + 1} = 'contains indirect cells'; %#ok<AGROW>
        end
        if ~isempty(hMainDirect)
            legendHandles(end + 1) = hMainDirect; %#ok<AGROW>
            legendLabels{end + 1} = sprintf('main-room cfos (%s mean)', aggregationMode); %#ok<AGROW>
        end
        if ~isempty(hControlRest)
            legendHandles(end + 1) = hControlRest; %#ok<AGROW>
            legendLabels{end + 1} = sprintf('control non-fos (%s mean)', aggregationMode); %#ok<AGROW>
        end
        if ~isempty(hControlDirect)
            legendHandles(end + 1) = hControlDirect; %#ok<AGROW>
            legendLabels{end + 1} = sprintf('control cfos (%s mean)', aggregationMode); %#ok<AGROW>
        end
        if ~isempty(legendHandles)
            legend(ax, legendHandles, legendLabels, 'Location', 'best');
        end
    end

    statsRows = table;
    statsRows = [statsRows; build_comparison_row(stateName, featureName, ...
        sprintf('main-room non-fos (%s mean)', aggregationMode), sprintf('main-room cfos (%s mean)', aggregationMode), ...
        mainRestValues, mainDirectValues, nnz(mainRestHasIndirect), pMain, zMain)]; %#ok<AGROW>
    statsRows = [statsRows; build_comparison_row(stateName, featureName, ...
        sprintf('control non-fos (%s mean)', aggregationMode), sprintf('control cfos (%s mean)', aggregationMode), ...
        controlRestValues, controlDirectValues, nnz(controlRestHasIndirect), pControl, zControl)]; %#ok<AGROW>
    statsRows = [statsRows; build_comparison_row(stateName, featureName, ...
        sprintf('main-room cfos (%s mean)', aggregationMode), sprintf('control cfos (%s mean)', aggregationMode), ...
        mainDirectValues, controlDirectValues, 0, pDirectRoom, zDirectRoom)]; %#ok<AGROW>
    statsRows = [statsRows; build_comparison_row(stateName, featureName, ...
        sprintf('main-room non-fos (%s mean)', aggregationMode), sprintf('control non-fos (%s mean)', aggregationMode), ...
        mainRestValues, controlRestValues, nnz(mainRestHasIndirect) + nnz(controlRestHasIndirect), pRestRoom, zRestRoom)]; %#ok<AGROW>
end

function [aggregationLabel, fileTag, subsetLabel] = aggregation_mode_labels(aggregationMode)

    if strcmpi(aggregationMode, 'animal')
        aggregationLabel = 'animal means';
        fileTag = 'AnimalMeans';
        subsetLabel = 'cfos_vs_non_fos_animal_mean';
    else
        aggregationLabel = 'session means';
        fileTag = 'SessionMeans';
        subsetLabel = 'cfos_vs_non_fos_session_mean';
    end
end

function txt = sanitize_threshold_text(value)

    txt = regexprep(sprintf('%.3f', value), '0+$', '');
    txt = regexprep(txt, '\.$', '');
    txt = strrep(txt, '.', 'p');
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

function statsRow = build_comparison_row(stateName, featureName, group1Label, group2Label, group1Values, group2Values, nIndirectHighlighted, pVal, zVal)

    statsRow = table( ...
        string(stateName), ...
        string(featureName), ...
        string(group1Label), ...
        string(group2Label), ...
        nnz(isfinite(group1Values)), ...
        nnz(isfinite(group2Values)), ...
        median(group1Values, 'omitnan'), ...
        median(group2Values, 'omitnan'), ...
        mean(group1Values, 'omitnan'), ...
        mean(group2Values, 'omitnan'), ...
        nIndirectHighlighted, ...
        pVal, ...
        zVal, ...
        'VariableNames', {'State','Feature','Group1Label','Group2Label','NGroup1','NGroup2', ...
        'MedianGroup1','MedianGroup2','MeanGroup1','MeanGroup2','NIndirectHighlighted', ...
        'RankSumPValue','RankSumZValue'});
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

function plot_group_summary(ax, x0, vals, colorVal)

    vals = vals(isfinite(vals));
    if isempty(vals)
        return
    end

    medVal = median(vals, 'omitnan');
    q = prctile(vals, [25 75]);

    plot(ax, [x0 x0], q, '-', 'Color', colorVal, 'LineWidth', 2);
    plot(ax, [x0 - 0.12, x0 + 0.12], [medVal medVal], '-', 'Color', colorVal, 'LineWidth', 3);
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

function folderName = four_condition_class_folder_name(className)

    if strcmpi(className, 'interneuron')
        folderName = 'Interneurons';
    elseif strcmpi(className, 'place cells')
        folderName = 'Place cells';
    elseif strcmpi(className, 'non-place principal')
        folderName = 'Principal non-place';
    else
        folderName = 'Principal cells';
    end
end

function write_column_csv(values, filename)

    values = values(:);

    [folder,~,~] = fileparts(filename);
    if ~exist(folder, 'dir')
        mkdir(folder);
    end

    fid = fopen(filename, 'w');
    if fid == -1
        return
    end
    cleaner = onCleanup(@() fclose(fid)); %#ok<NASGU>

    for i = 1:numel(values)
        if isfinite(values(i))
            fprintf(fid, '%.15g\n', values(i));
        else
            fprintf(fid, '\n');
        end
    end
end

function exportGroup = main_room_context_place_field_export_group(finalLabels)

    finalLabels = string(finalLabels);
    exportGroup = repmat("unclassified", size(finalLabels));
    exportGroup(finalLabels == "direct") = "main_room_cfos";
    exportGroup(finalLabels == "untagged" | finalLabels == "indirect") = "main_room_non_fos";
end

function bins = main_room_context_place_field_bin(fieldCounts)

    roundedCounts = round(fieldCounts(:));
    bins = repmat("unbinned", numel(roundedCounts), 1);
    bins(roundedCounts == 1) = "1 place field";
    bins(roundedCounts == 2) = "2 place fields";
    bins(roundedCounts >= 3) = "3+ place fields";
end

function [nValues, nCells, nOne, nTwo, nThreePlus, fOne, fTwo, fThreePlus] = ...
    count_main_room_context_place_field_bins(IncludedValues)

    fieldCounts = round(IncludedValues.PlaceFieldNumber);
    nValues = numel(fieldCounts);
    nCells = count_unique_main_room_context_cells(IncludedValues);
    nOne = nnz(fieldCounts == 1);
    nTwo = nnz(fieldCounts == 2);
    nThreePlus = nnz(fieldCounts >= 3);
    fOne = fraction_or_nan(nOne, nValues);
    fTwo = fraction_or_nan(nTwo, nValues);
    fThreePlus = fraction_or_nan(nThreePlus, nValues);
end

function [groupMeans, hasIndirect] = aggregate_group_means(T, values, includeMask, indirectMask, aggregationMode)

    groupMeans = [];
    hasIndirect = false(0, 1);

    aggregateIds = get_aggregation_ids(T, aggregationMode);
    validMask = includeMask & isfinite(values) & aggregateIds ~= "";
    if ~any(validMask)
        return
    end

    ids = aggregateIds(validMask);
    x = values(validMask);
    indirectIds = unique(aggregateIds(indirectMask & isfinite(values) & aggregateIds ~= ""), 'stable');

    [uniqueIds, ~, ic] = unique(ids, 'stable');
    groupMeans = nan(numel(uniqueIds), 1);
    hasIndirect = false(numel(uniqueIds), 1);

    for i = 1:numel(uniqueIds)
        memberMask = ic == i;
        groupMeans(i) = mean(x(memberMask), 'omitnan');
        hasIndirect(i) = any(indirectIds == uniqueIds(i));
    end
end

function [hMain, hIndirect] = scatter_aggregated_group(ax, xCenter, values, hasIndirect, faceColor)

    hMain = gobjects(0);
    hIndirect = gobjects(0);

    if isempty(values)
        return
    end

    hasIndirect = logical(hasIndirect(:));
    values = values(:);

    plainMask = ~hasIndirect;
    indirectMask = hasIndirect;

    if any(plainMask)
        x = xCenter + 0.16 * (rand(nnz(plainMask), 1) - 0.5);
        hMain = scatter(ax, x, values(plainMask), 32, faceColor, 'filled', ...
            'MarkerFaceAlpha', 0.80, 'MarkerEdgeColor', [0.1 0.1 0.1], 'LineWidth', 0.8);
    end

    if any(indirectMask)
        x = xCenter + 0.16 * (rand(nnz(indirectMask), 1) - 0.5);
        hIndirect = scatter(ax, x, values(indirectMask), 34, faceColor, 'filled', ...
            'MarkerFaceAlpha', 0.80, 'MarkerEdgeColor', [0.15 0.55 0.15], 'LineWidth', 1.0);
    end
end

function nCells = count_unique_main_room_context_cells(IncludedValues)

    if isempty(IncludedValues) || height(IncludedValues) == 0
        nCells = 0;
        return
    end

    cellKeys = strcat(string(IncludedValues.SessionIndex), "_", string(IncludedValues.Cell));
    nCells = numel(unique(cellKeys));
end

function value = fraction_or_nan(numerator, denominator)

    if denominator > 0
        value = numerator / denominator;
    else
        value = NaN;
    end
end

function aggregateIds = get_aggregation_ids(T, aggregationMode)

    if strcmpi(aggregationMode, 'animal')
        aggregateIds = string(T.AnimalID);
    else
        aggregateIds = string(T.SessionIndex);
    end
end
