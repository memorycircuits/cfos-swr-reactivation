function out = plot_interneuron_quantification_figures(workbookFile, varargin)
%PLOT_INTERNEURON_QUANTIFICATION_FIGURES Plot interneuron marker distributions.
%   OUT = PLOT_INTERNEURON_QUANTIFICATION_FIGURES(WORKBOOKFILE) imports
%   the interneuron quantification workbook and writes heatmaps, bar plots,
%   distribution profiles, summary CSV files, a slice/stack audit workbook,
%   an excluded-row audit workbook, and plotting notes.
%
%   WORKBOOKFILE is the Quantifications.xlsx file. If omitted, the function
%   searches the Results and SDC microscopy folders.
%
%   Name-value inputs:
%     OutputFolder             Folder for figures and exported tables.
%     SaveFigures              Write PNG/PDF/CSV/XLSX outputs when true.
%     CloseFigures             Close generated MATLAB figure handles.
%     IncludeNegativeMarkers   Include negative marker classes in plots.
%     NormalizeCountsBySlices  Normalize to 60x-equivalent stacks.
%     CacheFile                MAT cache path for the imported workbook.
%     ReimportExcel            Ignore the cache and re-read the workbook.
%
%   OUT contains the imported data, plotted summaries, slice audit tables,
%   figure index, shared bar-axis limits, output folder, and figure handles.

if nargin < 1
    workbookFile = [];
elseif is_name_value_start(workbookFile)
    varargin = [{workbookFile}, varargin];
    workbookFile = [];
end

parser = inputParser;
addParameter(parser, 'OutputFolder', default_output_folder(), @(x) ischar(x) || isstring(x));
addParameter(parser, 'SaveFigures', true, @(x) islogical(x) && isscalar(x));
addParameter(parser, 'CloseFigures', false, @(x) islogical(x) && isscalar(x));
addParameter(parser, 'IncludeNegativeMarkers', false, @(x) islogical(x) && isscalar(x));
addParameter(parser, 'NormalizeCountsBySlices', true, @(x) islogical(x) && isscalar(x));
addParameter(parser, 'CacheFile', '', @(x) ischar(x) || isstring(x));
addParameter(parser, 'ReimportExcel', false, @(x) islogical(x) && isscalar(x));
parse(parser, varargin{:});
opts = parser.Results;
opts.OutputFolder = char(opts.OutputFolder);

data = import_interneuron_quantifications(workbookFile, 'IncludeZeroCounts', true, ...
    'CacheFile', opts.CacheFile, 'ReimportExcel', opts.ReimportExcel);
warn_if_expected_denominators_missing(data);

animalRegion = filter_plot_markers(data.byAnimalRegion, opts.IncludeNegativeMarkers);
datasetRegion = filter_plot_markers(data.byDatasetMarkerRegion, opts.IncludeNegativeMarkers);
sliceAudit = make_slice_audit_tables(filter_plot_markers(data.byAnimalRegion, false));
barAxisLimits = make_bar_axis_limits(animalRegion, opts.NormalizeCountsBySlices);

if isempty(animalRegion)
    error('InterneuronPlot:NoData', 'No interneuron quantification rows were available for plotting.');
end

if opts.SaveFigures && ~exist(opts.OutputFolder, 'dir')
    mkdir(opts.OutputFolder);
end

datasets = unique(animalRegion(:, {'Dataset', 'MouseLine', 'Reporter'}), 'rows');
panels = unique(animalRegion(:, {'Dataset', 'MouseLine', 'Reporter', 'MarkerPanel'}), 'rows');
figures = gobjects(0, 1);
figureInfo = table(strings(0, 1), strings(0, 1), strings(0, 1), strings(0, 1), ...
    'VariableNames', {'Dataset', 'MarkerPanel', 'FigureType', 'FileBase'});

for datasetIdx = 1:height(datasets)
    datasetMask = animalRegion.Dataset == datasets.Dataset(datasetIdx);
    datasetAnimalRegion = animalRegion(datasetMask, :);
    if isempty(datasetAnimalRegion)
        continue
    end

    filePrefix = sanitize_filename(sprintf('%s_all_markers', ...
        char(datasets.Dataset(datasetIdx))));
    fig = plot_answer_heatmap(datasetAnimalRegion, datasets(datasetIdx, :), ...
        opts.NormalizeCountsBySlices);
    figures(end + 1, 1) = fig; %#ok<AGROW>
    figureInfo = add_figure_info(figureInfo, datasets(datasetIdx, :), ...
        "all_marker_panels", "answer_heatmap", filePrefix + "_answer_heatmap");

    selectedDatasetAnimalRegion = filter_selected_heatmap_markers(datasetAnimalRegion);
    if ~isempty(selectedDatasetAnimalRegion)
        selectedFilePrefix = sanitize_filename(sprintf('%s_selected_markers', ...
            char(datasets.Dataset(datasetIdx))));
        fig = plot_answer_heatmap(selectedDatasetAnimalRegion, datasets(datasetIdx, :), ...
            opts.NormalizeCountsBySlices, 1.2, "selected", ...
            "selected PV/SOM/proCCK/VIP classes");
        figures(end + 1, 1) = fig; %#ok<AGROW>
        figureInfo = add_figure_info(figureInfo, datasets(datasetIdx, :), ...
            "PV/SOM + proCCK/VIP selected", "selected_answer_heatmap", ...
            selectedFilePrefix + "_answer_heatmap");
    end

    fig = plot_total_distribution_profiles(datasetAnimalRegion, datasets(datasetIdx, :), ...
        opts.NormalizeCountsBySlices);
    figures(end + 1, 1) = fig; %#ok<AGROW>
    figureInfo = add_figure_info(figureInfo, datasets(datasetIdx, :), ...
        "all_marker_panels", "total_distribution_profiles", ...
        filePrefix + "_total_distribution_profiles");
end

for panelIdx = 1:height(panels)
    panelMask = animalRegion.Dataset == panels.Dataset(panelIdx) & ...
        animalRegion.MarkerPanel == panels.MarkerPanel(panelIdx);

    panelAnimalRegion = animalRegion(panelMask, :);

    if isempty(panelAnimalRegion)
        continue
    end

    filePrefix = sanitize_filename(sprintf('%s_%s', ...
        char(panels.Dataset(panelIdx)), char(panels.MarkerPanel(panelIdx))));

    xAxisMax = get_bar_axis_max(barAxisLimits, panels.Dataset(panelIdx));
    fig = plot_total_active_bars(panelAnimalRegion, panels(panelIdx, :), ...
        opts.NormalizeCountsBySlices, xAxisMax);
    figures(end + 1, 1) = fig; %#ok<AGROW>
    figureInfo = add_figure_info(figureInfo, panels(panelIdx, :), ...
        panels.MarkerPanel(panelIdx), "total_active_bars", filePrefix + "_total_active_bars");
end

if opts.SaveFigures
    writetable(animalRegion, fullfile(opts.OutputFolder, 'interneuron_by_animal_region.csv'));
    writetable(data.byAnimalMarker, fullfile(opts.OutputFolder, 'interneuron_by_animal_marker.csv'));
    writetable(datasetRegion, fullfile(opts.OutputFolder, 'interneuron_by_dataset_marker_region.csv'));
    writetable(data.totalAvailability, fullfile(opts.OutputFolder, 'interneuron_denominator_availability.csv'));
    writetable(figureInfo, fullfile(opts.OutputFolder, 'interneuron_figure_index.csv'));
    writetable(barAxisLimits, fullfile(opts.OutputFolder, 'interneuron_bar_axis_limits.csv'));
    write_slice_audit_workbook(opts.OutputFolder, sliceAudit);
    write_excluded_rows_workbook(opts.OutputFolder, data);

    for figIdx = 1:numel(figures)
        save_figure(figures(figIdx), opts.OutputFolder, char(figureInfo.FileBase(figIdx)));
    end
    write_plot_notes(opts.OutputFolder, data, figureInfo, opts, barAxisLimits);
    fprintf('Interneuron quantification outputs saved to:\n  %s\n', opts.OutputFolder);
    fprintf('Slice/stack audit workbook:\n  %s\n', ...
        fullfile(opts.OutputFolder, 'interneuron_analyzed_slices.xlsx'));
    fprintf('Excluded-row audit workbook:\n  %s\n', ...
        fullfile(opts.OutputFolder, 'interneuron_excluded_rows.xlsx'));
elseif ~opts.SaveFigures
    fprintf('SaveFigures is false. No files were written. OutputFolder is:\n  %s\n', ...
        opts.OutputFolder);
end

if opts.CloseFigures
    close(figures(ishandle(figures)));
end

out = struct();
out.data = data;
out.animalRegion = animalRegion;
out.datasetRegion = datasetRegion;
out.sliceAudit = sliceAudit;
out.excludedSlices = data.excludedSlices;
out.excludedRows = data.excludedRows;
out.barAxisLimits = barAxisLimits;
out.figureInfo = figureInfo;
out.figures = figures;
out.outputFolder = string(opts.OutputFolder);
out.normalizeCountsBySlices = opts.NormalizeCountsBySlices;
end

function T = filter_plot_markers(T, includeNegativeMarkers)
if includeNegativeMarkers || isempty(T)
    return
end
T = T(~contains(T.Marker, "-"), :);
end

function T = filter_selected_heatmap_markers(T)
if isempty(T)
    return
end
pvSomMarkers = ["PV", "SST/SOM", "PV+SST/SOM"];
proCckVipMarkers = ["proCCK", "VIP"];
keepRows = (T.MarkerPanel == "PV/SOM" & ismember(T.Marker, pvSomMarkers)) | ...
    (T.MarkerPanel == "proCCK/VIP" & ismember(T.Marker, proCckVipMarkers));
T = T(keepRows, :);
end

function warn_if_expected_denominators_missing(data)
if ~isfield(data, 'totalAvailability') || isempty(data.totalAvailability)
    return
end

availability = data.totalAvailability;
expectedRows = (availability.Dataset == "cfos_TRE_mKate" & ...
    (availability.MarkerPanel == "PV/SOM" | availability.MarkerPanel == "proCCK/VIP")) | ...
    (availability.Dataset == "intrinsic_cFos" & ...
    (availability.MarkerPanel == "PV/SOM" | availability.MarkerPanel == "SNCG/CCK8" | ...
    availability.MarkerPanel == "proCCK/VIP"));

missingRows = expectedRows & availability.NRows > 0 & ...
    availability.RowsWithRegionDenominator == 0 & ...
    availability.RowsWithAllRegionDenominator == 0 & ...
    availability.RowsWithMissingDenominator == availability.NRows;

if any(missingRows)
    badRows = availability(missingRows, :);
    labels = strings(height(badRows), 1);
    for rowIdx = 1:height(badRows)
        labels(rowIdx) = badRows.Dataset(rowIdx) + " / " + ...
            badRows.MarkerPanel(rowIdx) + " / " + badRows.Marker(rowIdx);
    end
    warning('InterneuronPlot:MissingExpectedDenominators', ...
        ['Some marker groups have no imported denominators: %s. ', ...
        'Check that the current importer is on the MATLAB path, clear old functions, ', ...
        'and run with ReimportExcel=true or delete the stale cache file.'], ...
        char(strjoin(labels, '; ')));
end
end

function fig = plot_answer_heatmap(T, datasetInfo, normalizeCountsBySlices, widthScale, markerLabelMode, titleSuffix)
if nargin < 4 || isempty(widthScale)
    widthScale = 1;
end
if nargin < 5 || strlength(markerLabelMode) == 0
    markerLabelMode = "source";
end
if nargin < 6
    titleSuffix = "";
end

regions = ordered_regions(T.Region);
rowInfo = make_combined_marker_rows(T);
[totalMatrix, activeMatrix, fractionMatrix, labelMatrix, fallbackMatrix] = ...
    make_answer_matrices(T, rowInfo, regions, normalizeCountsBySlices);

countUnit = count_unit_label(normalizeCountsBySlices);
countFormat = count_value_format(normalizeCountsBySlices);

markerLabels = heatmap_marker_labels(rowInfo);
markerLabels = display_heatmap_marker_labels(markerLabels, markerLabelMode);
markerLabels(any(fallbackMatrix, 2)) = markerLabels(any(fallbackMatrix, 2)) + " *";
regionLabels = anatomical_region_labels(regions);
figWidth = max(22, min(30, 18 + height(rowInfo) * 0.75)) * widthScale;
figHeight = max(29, min(36, 10 + numel(regions) * 5.0));
fig = figure('Color', 'w', 'Units', 'centimeters', 'Position', [2 2 figWidth figHeight]);
tiledlayout(fig, 3, 1, 'Padding', 'compact', 'TileSpacing', 'compact');

ax = nexttile;
plot_numeric_heatmap(ax, totalMatrix', regionLabels, markerLabels, ...
    'where are the cell types?', char("mean total " + countUnit), countFormat, widthScale);

ax = nexttile;
plot_numeric_heatmap(ax, activeMatrix', regionLabels, markerLabels, ...
    sprintf('how many are %s+?', char(datasetInfo.Reporter)), ...
    sprintf('mean %s+ %s', char(datasetInfo.Reporter), char(countUnit)), countFormat, widthScale);

ax = nexttile;
plot_numeric_heatmap(ax, fractionMatrix', regionLabels, markerLabels, ...
    sprintf('what fraction are %s+?', char(datasetInfo.Reporter)), ...
    sprintf('mean %% %s+', char(datasetInfo.Reporter)), '%.0f%%', widthScale);
overlay_fallback_marks(ax, fallbackMatrix');

titleText = string(sprintf('%s | %s: location, active count, active fraction', ...
    char(datasetInfo.Dataset), char(datasetInfo.Reporter)));
if strlength(titleSuffix) > 0
    titleText = titleText + " | " + string(titleSuffix);
end
sgtitle(fig, char(titleText), 'Interpreter', 'none');

if any(fallbackMatrix(:))
    annotation(fig, 'textbox', [0.02 0.01 0.96 0.04], ...
        'String', ['* no region-level total in workbook: total panel is blank; ', ...
        'fraction uses active/all-region marker total where available'], ...
        'EdgeColor', 'none', 'FontSize', 8, 'Interpreter', 'none', ...
        'Color', [0.35 0.05 0.05]);
end
end

function labels = display_heatmap_marker_labels(labels, markerLabelMode)
if markerLabelMode ~= "selected"
    return
end
labels(labels == "SST/SOM") = "SOM";
labels(labels == "PV+SST/SOM") = "PV/SOM";
end

function [totalMatrix, activeMatrix, fractionMatrix, labelMatrix, fallbackMatrix] = ...
    make_answer_matrices(T, rowInfo, regions, normalizeCountsBySlices)
countVars = count_variable_names(normalizeCountsBySlices);
totalMatrix = nan(height(rowInfo), numel(regions));
activeMatrix = nan(height(rowInfo), numel(regions));
fractionMatrix = nan(height(rowInfo), numel(regions));
labelMatrix = strings(height(rowInfo), numel(regions));
fallbackMatrix = false(height(rowInfo), numel(regions));

for rowIdx = 1:height(rowInfo)
    for regionIdx = 1:numel(regions)
        rowMask = T.MarkerPanel == rowInfo.MarkerPanel(rowIdx) & ...
            T.Marker == rowInfo.Marker(rowIdx) & T.Region == regions(regionIdx);
        rows = T(rowMask, :);
        if isempty(rows)
            continue
        end

        activeMatrix(rowIdx, regionIdx) = mean_finite(rows.(countVars.Active));
        regionTotal = mean_finite(rows.(countVars.TotalRegion));
        allRegionTotal = mean_finite(rows.(countVars.TotalAllRegions));

        if isfinite(regionTotal)
            totalMatrix(rowIdx, regionIdx) = regionTotal;
            fractionMatrix(rowIdx, regionIdx) = mean_finite(rows.ActivePercentWithinRegion);
            labelMatrix(rowIdx, regionIdx) = format_count_pair( ...
                activeMatrix(rowIdx, regionIdx), totalMatrix(rowIdx, regionIdx), "");
        else
            fallbackMatrix(rowIdx, regionIdx) = true;
            fractionMatrix(rowIdx, regionIdx) = mean_finite(rows.ActivePercentOfAllRegions);
            if isfinite(allRegionTotal)
                labelMatrix(rowIdx, regionIdx) = format_count_pair( ...
                    activeMatrix(rowIdx, regionIdx), allRegionTotal, "*");
            else
                labelMatrix(rowIdx, regionIdx) = format_count_value(activeMatrix(rowIdx, regionIdx)) + "/?";
            end
        end
    end
end
end

function labels = heatmap_marker_labels(rowInfo)
labels = rowInfo.Marker;
if numel(unique(labels)) < numel(labels)
    labels = rowInfo.RowLabel;
end
end

function plot_numeric_heatmap(ax, values, rowLabels, colLabels, titleText, colorbarText, labelFormat, widthScale)
if nargin < 8 || isempty(widthScale)
    widthScale = 1;
end
imageHandle = imagesc(ax, values);
set(imageHandle, 'AlphaData', ~isnan(values));
finiteValues = values(isfinite(values));
if isempty(finiteValues)
    colorMax = 1;
else
    colorMax = max(finiteValues);
    if colorMax == 0
        colorMax = 1;
    end
end
set(ax, 'Color', [0.92 0.92 0.92], 'CLim', [0 colorMax], ...
    'XTick', 1:numel(colLabels), 'XTickLabel', cellstr(colLabels), ...
    'YTick', 1:numel(rowLabels), 'YTickLabel', cellstr(rowLabels), ...
    'YDir', 'reverse', 'TickDir', 'out');
xtickangle(ax, 45);
pbaspect(ax, [widthScale 1 1]);
colormap(ax, parula);
colorHandle = colorbar(ax);
ylabel(colorHandle, colorbarText, 'Interpreter', 'none');
title(ax, titleText, 'Interpreter', 'none');
draw_heatmap_values(ax, values, labelFormat);
end

function draw_heatmap_values(ax, values, labelFormat)
for rowIdx = 1:size(values, 1)
    for colIdx = 1:size(values, 2)
        value = values(rowIdx, colIdx);
        if ~isfinite(value)
            text(ax, colIdx, rowIdx, 'n/a', 'HorizontalAlignment', 'center', ...
                'FontSize', 7, 'Color', [0.35 0.35 0.35], 'Interpreter', 'none');
            continue
        end
        if contains(labelFormat, '%%')
            labelText = sprintf(labelFormat, value);
        else
            labelText = sprintf(labelFormat, value);
        end
        textColor = [0.08 0.08 0.08];
        clim = get(ax, 'CLim');
        if value > mean(clim)
            textColor = [1 1 1];
        end
        text(ax, colIdx, rowIdx, labelText, 'HorizontalAlignment', 'center', ...
            'FontSize', 7, 'Color', textColor, 'Interpreter', 'none');
    end
end
end

function overlay_fallback_marks(ax, fallbackMatrix)
for rowIdx = 1:size(fallbackMatrix, 1)
    for colIdx = 1:size(fallbackMatrix, 2)
        if fallbackMatrix(rowIdx, colIdx)
            text(ax, colIdx + 0.35, rowIdx - 0.32, '*', ...
                'HorizontalAlignment', 'center', 'FontSize', 9, ...
                'FontWeight', 'bold', 'Color', [0.80 0.02 0.02]);
        end
    end
end
end

function fig = plot_total_distribution_profiles(T, datasetInfo, normalizeCountsBySlices)
regions = ordered_regions(T.Region);
rowInfo = make_combined_marker_rows(T);
markerLabels = heatmap_marker_labels(rowInfo);
colors = marker_colors(markerLabels);
countUnit = count_unit_label(normalizeCountsBySlices);
y = 1:numel(regions);
xAxisMax = 1;

fig = figure('Color', 'w', 'Units', 'centimeters', 'Position', [2 2 23 14]);
ax = axes(fig);
hold(ax, 'on');

legendHandles = gobjects(0, 1);
legendLabels = strings(0, 1);
for markerIdx = 1:height(rowInfo)
    [meanCounts, semCounts, nMice] = marker_distribution_stats(T, rowInfo(markerIdx, :), ...
        regions, normalizeCountsBySlices);
    finiteRows = isfinite(meanCounts);
    if ~any(finiteRows)
        continue
    end

    color = colors(markerIdx, :);
    lineHandle = plot(ax, meanCounts(finiteRows), y(finiteRows), '-o', ...
        'Color', color, 'MarkerFaceColor', color, 'MarkerEdgeColor', [0.05 0.05 0.05], ...
        'LineWidth', 1.3, 'MarkerSize', 5);
    draw_horizontal_errorbars(ax, meanCounts, y, semCounts, color);
    finiteExtent = meanCounts + semCounts;
    finiteExtent = finiteExtent(isfinite(finiteExtent));
    if ~isempty(finiteExtent)
        xAxisMax = max(xAxisMax, max(finiteExtent));
    end

    legendHandles(end + 1, 1) = lineHandle; %#ok<AGROW>
    legendLabels(end + 1, 1) = markerLabels(markerIdx) + " (" + mouse_count_label(nMice) + ")"; %#ok<AGROW>
end

xlim(ax, [0 nice_axis_limit(xAxisMax)]);
set(ax, 'YLim', [0.5 numel(regions) + 0.5], 'YTick', y, ...
    'YTickLabel', cellstr(anatomical_region_labels(regions)), ...
    'YDir', 'reverse', 'TickDir', 'out');
xlabel(ax, char("mean total " + countUnit));
ylabel(ax, 'Hippocampal layer');
title(ax, sprintf('%s | %s: total interneuron distribution across layers', ...
    char(datasetInfo.Dataset), char(datasetInfo.Reporter)), 'Interpreter', 'none');
if ~isempty(legendHandles)
    legend(ax, legendHandles, cellstr(legendLabels), 'Location', 'eastoutside', ...
        'Interpreter', 'none');
end
box(ax, 'off');
hold(ax, 'off');
end

function [meanCounts, semCounts, nMice] = marker_distribution_stats(T, markerRow, regions, normalizeCountsBySlices)
countVars = count_variable_names(normalizeCountsBySlices);
markerMask = T.MarkerPanel == markerRow.MarkerPanel & T.Marker == markerRow.Marker;
mice = unique(T.Mouse(markerMask & isfinite(T.Mouse)));
profiles = nan(numel(mice), numel(regions));

for mouseIdx = 1:numel(mice)
    regionValues = nan(1, numel(regions));
    for regionIdx = 1:numel(regions)
        rowMask = markerMask & T.Mouse == mice(mouseIdx) & T.Region == regions(regionIdx);
        rows = T(rowMask, :);
        if isempty(rows)
            continue
        end
        regionValues(regionIdx) = mean_finite(rows.(countVars.TotalRegion));
    end

    profiles(mouseIdx, :) = regionValues;
end

validProfiles = any(isfinite(profiles), 2);
nMice = sum(validProfiles);
meanCounts = nan(1, numel(regions));
semCounts = nan(1, numel(regions));
for regionIdx = 1:numel(regions)
    values = profiles(:, regionIdx);
    values = values(isfinite(values));
    if isempty(values)
        continue
    end
    meanCounts(regionIdx) = mean(values);
    semCounts(regionIdx) = sem_finite(values);
end
end

function draw_horizontal_errorbars(ax, xValues, yValues, errors, color)
capHalfHeight = 0.06;
for idx = 1:numel(xValues)
    if ~isfinite(xValues(idx)) || ~isfinite(errors(idx))
        continue
    end
    xLow = max(0, xValues(idx) - errors(idx));
    xHigh = xValues(idx) + errors(idx);
    line(ax, [xLow xHigh], [yValues(idx) yValues(idx)], ...
        'Color', color, 'LineWidth', 1.0);
    line(ax, [xLow xLow], [yValues(idx) - capHalfHeight yValues(idx) + capHalfHeight], ...
        'Color', color, 'LineWidth', 1.0);
    line(ax, [xHigh xHigh], [yValues(idx) - capHalfHeight yValues(idx) + capHalfHeight], ...
        'Color', color, 'LineWidth', 1.0);
end
end

function value = sem_finite(values)
values = values(isfinite(values));
if numel(values) <= 1
    value = 0;
else
    value = std(values) / sqrt(numel(values));
end
end

function fig = plot_total_active_bars(T, panel, normalizeCountsBySlices, xAxisMax)
if nargin < 4 || ~isfinite(xAxisMax) || xAxisMax <= 0
    xAxisMax = compute_bar_axis_max(T, normalizeCountsBySlices);
end
regions = ordered_regions(T.Region);
markers = ordered_markers(T.Marker, panel.MarkerPanel);
colors = marker_colors(markers);
countVars = count_variable_names(normalizeCountsBySlices);
countUnit = count_unit_label(normalizeCountsBySlices);
nMarkers = numel(markers);
nCols = min(3, nMarkers);
nRows = ceil(nMarkers / nCols);

fig = figure('Color', 'w', 'Units', 'centimeters', ...
    'Position', [2 2 8.0 * nCols, 5.9 * nRows + 1.6]);
tiledlayout(fig, nRows, nCols, 'Padding', 'compact', 'TileSpacing', 'compact');

for markerIdx = 1:nMarkers
    ax = nexttile;
    hold(ax, 'on');
    totalValues = nan(1, numel(regions));
    activeValues = nan(1, numel(regions));
    fallbackValues = false(1, numel(regions));
    markerRows = T(T.Marker == markers(markerIdx), :);
    nMice = count_unique_finite(markerRows.Mouse);
    sliceCounts = mouse_slice_counts(markerRows);
    fovFactors = mouse_fov_area_factors(markerRows);
    normalizationCounts = mouse_normalization_factors(markerRows);

    for regionIdx = 1:numel(regions)
        rowMask = T.Marker == markers(markerIdx) & T.Region == regions(regionIdx);
        rows = T(rowMask, :);
        activeValues(regionIdx) = mean_finite(rows.(countVars.Active));
        totalValues(regionIdx) = mean_finite(rows.(countVars.TotalRegion));
        if ~isfinite(totalValues(regionIdx))
            fallbackValues(regionIdx) = true;
        end
    end

    y = 1:numel(regions);
    lightColor = lighten_color(colors(markerIdx, :), 0.62);
    totalBarValues = totalValues;
    totalBarValues(~isfinite(totalBarValues)) = 0;
    activeBarValues = activeValues;
    activeBarValues(~isfinite(activeBarValues)) = 0;

    barh(ax, y, totalBarValues, 0.68, 'FaceColor', lightColor, ...
        'EdgeColor', [0.25 0.25 0.25], 'LineWidth', 0.45);
    barh(ax, y, activeBarValues, 0.42, 'FaceColor', colors(markerIdx, :), ...
        'EdgeColor', [0.05 0.05 0.05], 'LineWidth', 0.35);

    for regionIdx = 1:numel(regions)
        xRight = max(totalBarValues(regionIdx), activeBarValues(regionIdx));
        if xRight == 0
            xRight = xAxisMax * 0.02;
        end
        [labelX, horizontalAlignment] = bar_label_position(xRight, xAxisMax);
        if fallbackValues(regionIdx)
            text(ax, labelX, y(regionIdx), ...
                char(format_count_value(activeBarValues(regionIdx)) + "/?"), ...
                'HorizontalAlignment', horizontalAlignment, 'VerticalAlignment', 'middle', ...
                'FontSize', 7, 'Color', [0.65 0 0]);
        else
            text(ax, labelX, y(regionIdx), ...
                char(format_count_pair(activeBarValues(regionIdx), totalBarValues(regionIdx), "")), ...
                'HorizontalAlignment', horizontalAlignment, 'VerticalAlignment', 'middle', ...
                'FontSize', 7, 'Color', [0.1 0.1 0.1]);
        end
    end

    xlim(ax, [0 xAxisMax]);
    set(ax, 'YLim', [0.5 numel(regions) + 0.5], 'YTick', y, ...
        'YTickLabel', cellstr(anatomical_region_labels(regions)), ...
        'YDir', 'reverse', 'TickDir', 'out');
    xLabelLines = {char("mean " + countUnit)};
    if normalizeCountsBySlices
        xLabelLines{end + 1} = char("stacks/mouse: " + format_slice_summary(sliceCounts));
        if any(fovFactors(isfinite(fovFactors)) > 1)
            xLabelLines{end + 1} = char("FOV area factor: " + ...
                format_factor_summary(fovFactors) + "; 60x-eq stacks/mouse: " + ...
                format_slice_summary(normalizationCounts));
        end
    end
    xlabel(ax, xLabelLines, 'FontSize', 8);
    title(ax, char(markers(markerIdx) + " (" + mouse_count_label(nMice) + ")"), ...
        'Interpreter', 'none');
    box(ax, 'off');
    hold(ax, 'off');
end

sgtitle(fig, sprintf('%s | %s | %s: total cells with active/tagged overlay', ...
    char(panel.Dataset), char(panel.MarkerPanel), char(panel.Reporter)), ...
    'Interpreter', 'none');
end

function rowInfo = make_combined_marker_rows(T)
panelOrder = ordered_marker_panels(T.MarkerPanel);
MarkerPanel = strings(0, 1);
Marker = strings(0, 1);
RowLabel = strings(0, 1);

for panelIdx = 1:numel(panelOrder)
    panel = panelOrder(panelIdx);
    panelRows = T.MarkerPanel == panel;
    markers = ordered_markers(T.Marker(panelRows), panel);
    for markerIdx = 1:numel(markers)
        MarkerPanel(end + 1, 1) = panel; %#ok<AGROW>
        Marker(end + 1, 1) = markers(markerIdx); %#ok<AGROW>
        RowLabel(end + 1, 1) = panel + " | " + markers(markerIdx); %#ok<AGROW>
    end
end

rowInfo = table(MarkerPanel, Marker, RowLabel);
end

function panels = ordered_marker_panels(panelValues)
preferred = ["PV/SOM"; "SNCG/CCK8"; "proCCK/VIP"];
panels = strings(0, 1);
for idx = 1:numel(preferred)
    if any(panelValues == preferred(idx))
        panels(end + 1, 1) = preferred(idx); %#ok<AGROW>
    end
end
remaining = unique(panelValues);
remaining = remaining(~ismember(remaining, panels));
panels = [panels; remaining(:)];
end

function value = get_region_mean_value(T, marker, region, valueName)
rowMask = T.Marker == marker & T.Region == region;
if ~any(rowMask)
    value = nan;
    return
end
values = T.(valueName)(rowMask);
value = mean_finite(values);
end

function value = mean_finite(values)
values = values(isfinite(values));
if isempty(values)
    value = nan;
else
    value = mean(values);
end
end

function axisLimits = make_bar_axis_limits(T, normalizeCountsBySlices)
datasets = unique(T(:, {'Dataset', 'MouseLine', 'Reporter'}), 'rows');
BarXAxisMax = nan(height(datasets), 1);

for datasetIdx = 1:height(datasets)
    rows = T(T.Dataset == datasets.Dataset(datasetIdx), :);
    BarXAxisMax(datasetIdx) = compute_bar_axis_max(rows, normalizeCountsBySlices);
end

axisLimits = [datasets, table(BarXAxisMax)];
end

function xAxisMax = get_bar_axis_max(axisLimits, datasetName)
rowMask = axisLimits.Dataset == datasetName;
if any(rowMask)
    xAxisMax = axisLimits.BarXAxisMax(find(rowMask, 1, 'first'));
else
    xAxisMax = nan;
end
end

function xAxisMax = compute_bar_axis_max(T, normalizeCountsBySlices)
countVars = count_variable_names(normalizeCountsBySlices);
regions = ordered_regions(T.Region);
panelOrder = ordered_marker_panels(T.MarkerPanel);
dataMax = 0;

for panelIdx = 1:numel(panelOrder)
    panel = panelOrder(panelIdx);
    panelRows = T.MarkerPanel == panel;
    markers = ordered_markers(T.Marker(panelRows), panel);
    for markerIdx = 1:numel(markers)
        for regionIdx = 1:numel(regions)
            rowMask = T.MarkerPanel == panel & T.Marker == markers(markerIdx) & ...
                T.Region == regions(regionIdx);
            rows = T(rowMask, :);
            if isempty(rows)
                continue
            end
            values = [mean_finite(rows.(countVars.Active)), ...
                mean_finite(rows.(countVars.TotalRegion))];
            values = values(isfinite(values));
            if ~isempty(values)
                dataMax = max(dataMax, max(values));
            end
        end
    end
end

if ~isfinite(dataMax) || dataMax <= 0
    dataMax = 1;
end
xAxisMax = nice_axis_limit(dataMax);
end

function limitValue = nice_axis_limit(value)
if ~isfinite(value) || value <= 0
    limitValue = 1;
    return
end

magnitude = 10 ^ floor(log10(value));
scaled = value / magnitude;
niceSteps = [1 1.25 1.5 2 2.5 3 4 5 7.5 10];
niceScaled = niceSteps(find(scaled <= niceSteps, 1, 'first'));
limitValue = niceScaled * magnitude;
end

function [xValue, horizontalAlignment] = bar_label_position(xRight, xAxisMax)
offset = xAxisMax * 0.025;
xValue = xRight + offset;
horizontalAlignment = 'left';
if xValue > xAxisMax * 0.98
    xValue = xAxisMax * 0.98;
    horizontalAlignment = 'right';
end
end

function audit = make_slice_audit_tables(T)
T = ensure_slice_columns(T);
audit = struct();
audit.datasetPanel = summarize_slice_groups(T, ...
    {'Dataset', 'MouseLine', 'Reporter', 'Sheet', 'MarkerPanel'});
audit.datasetMarker = summarize_slice_groups(T, ...
    {'Dataset', 'MouseLine', 'Reporter', 'Sheet', 'MarkerPanel', 'Marker'});
audit.mousePanel = summarize_slice_groups(T, ...
    {'Dataset', 'MouseLine', 'Reporter', 'Sheet', 'MarkerPanel', 'Mouse'});
audit.mouseMarkerRegion = mouse_marker_region_slice_table(T);
end

function T = ensure_slice_columns(T)
if isempty(T)
    return
end

if ~any(strcmp(T.Properties.VariableNames, 'NAnalyzedSlices'))
    if any(strcmp(T.Properties.VariableNames, 'NInputRows'))
        T.NAnalyzedSlices = T.NInputRows;
    else
        T.NAnalyzedSlices = nan(height(T), 1);
    end
end

if ~any(strcmp(T.Properties.VariableNames, 'NStacks'))
    T.NStacks = nan(height(T), 1);
end

if ~any(strcmp(T.Properties.VariableNames, 'MeanFovAreaFactor'))
    T.MeanFovAreaFactor = ones(height(T), 1);
end

if ~any(strcmp(T.Properties.VariableNames, 'CountNormalizationFactor'))
    T.CountNormalizationFactor = T.NAnalyzedSlices .* T.MeanFovAreaFactor;
end

if ~any(strcmp(T.Properties.VariableNames, 'NInputRows'))
    T.NInputRows = nan(height(T), 1);
end
end

function summary = summarize_slice_groups(T, groupVars)
summary = empty_slice_summary_table(groupVars);
if isempty(T)
    return
end

[G, groupTable] = findgroups(T(:, groupVars));
NMouse = nan(height(groupTable), 1);
MeanSlicesStacksPerMouse = nan(height(groupTable), 1);
MedianSlicesStacksPerMouse = nan(height(groupTable), 1);
MinSlicesStacksPerMouse = nan(height(groupTable), 1);
MaxSlicesStacksPerMouse = nan(height(groupTable), 1);
TotalMouseSlicesStacks = nan(height(groupTable), 1);
MeanFovAreaFactor = nan(height(groupTable), 1);
MeanCountNormalizationFactorPerMouse = nan(height(groupTable), 1);
NRows = nan(height(groupTable), 1);
MouseSliceStackCounts = strings(height(groupTable), 1);
MouseFovAreaFactors = strings(height(groupTable), 1);
MouseCountNormalizationFactors = strings(height(groupTable), 1);

for groupIdx = 1:height(groupTable)
    rows = T(G == groupIdx, :);
    sliceCounts = mouse_slice_counts(rows);
    fovFactors = mouse_fov_area_factors(rows);
    normalizationFactors = mouse_normalization_factors(rows);
    sliceCounts = sliceCounts(isfinite(sliceCounts));
    fovFactors = fovFactors(isfinite(fovFactors));
    normalizationFactors = normalizationFactors(isfinite(normalizationFactors));
    NMouse(groupIdx) = count_unique_finite(rows.Mouse);
    NRows(groupIdx) = height(rows);
    MouseSliceStackCounts(groupIdx) = format_mouse_slice_counts(rows);
    MouseFovAreaFactors(groupIdx) = format_mouse_value_counts(rows, 'MeanFovAreaFactor');
    MouseCountNormalizationFactors(groupIdx) = format_mouse_value_counts(rows, 'CountNormalizationFactor');
    if isempty(sliceCounts)
        continue
    end
    MeanSlicesStacksPerMouse(groupIdx) = mean(sliceCounts);
    MedianSlicesStacksPerMouse(groupIdx) = median(sliceCounts);
    MinSlicesStacksPerMouse(groupIdx) = min(sliceCounts);
    MaxSlicesStacksPerMouse(groupIdx) = max(sliceCounts);
    TotalMouseSlicesStacks(groupIdx) = sum(sliceCounts);
    if ~isempty(fovFactors)
        MeanFovAreaFactor(groupIdx) = mean(fovFactors);
    end
    if ~isempty(normalizationFactors)
        MeanCountNormalizationFactorPerMouse(groupIdx) = mean(normalizationFactors);
    end
end

summary = [groupTable, table(NMouse, MeanSlicesStacksPerMouse, ...
    MedianSlicesStacksPerMouse, MinSlicesStacksPerMouse, ...
    MaxSlicesStacksPerMouse, TotalMouseSlicesStacks, MeanFovAreaFactor, ...
    MeanCountNormalizationFactorPerMouse, MouseSliceStackCounts, ...
    MouseFovAreaFactors, MouseCountNormalizationFactors, NRows)];
end

function summary = empty_slice_summary_table(groupVars)
summary = table();
for varIdx = 1:numel(groupVars)
    varName = groupVars{varIdx};
    if strcmp(varName, 'Mouse')
        summary.(varName) = nan(0, 1);
    else
        summary.(varName) = strings(0, 1);
    end
end
summary.NMouse = nan(0, 1);
summary.MeanSlicesStacksPerMouse = nan(0, 1);
summary.MedianSlicesStacksPerMouse = nan(0, 1);
summary.MinSlicesStacksPerMouse = nan(0, 1);
summary.MaxSlicesStacksPerMouse = nan(0, 1);
summary.TotalMouseSlicesStacks = nan(0, 1);
summary.MeanFovAreaFactor = nan(0, 1);
summary.MeanCountNormalizationFactorPerMouse = nan(0, 1);
summary.MouseSliceStackCounts = strings(0, 1);
summary.MouseFovAreaFactors = strings(0, 1);
summary.MouseCountNormalizationFactors = strings(0, 1);
summary.NRows = nan(0, 1);
end

function T = mouse_marker_region_slice_table(T)
if isempty(T)
    T = table();
    return
end

keepVars = {'Dataset', 'MouseLine', 'Reporter', 'Sheet', 'MarkerPanel', ...
    'Mouse', 'Marker', 'Region', 'NAnalyzedSlices', 'MeanFovAreaFactor', ...
    'CountNormalizationFactor', 'NStacks', 'NInputRows'};
hasVar = ismember(keepVars, T.Properties.VariableNames);
T = unique(T(:, keepVars(hasVar)), 'rows');
end

function write_slice_audit_workbook(outputFolder, audit)
workbookFile = fullfile(outputFolder, 'interneuron_analyzed_slices.xlsx');
try
    if isfile(workbookFile)
        delete(workbookFile);
    end
    writetable(audit.datasetPanel, workbookFile, 'Sheet', 'DatasetPanelSummary');
    writetable(audit.datasetMarker, workbookFile, 'Sheet', 'MarkerSummary');
    writetable(audit.mousePanel, workbookFile, 'Sheet', 'MousePanelSlices');
    writetable(audit.mouseMarkerRegion, workbookFile, 'Sheet', 'MouseMarkerRegion');
catch err
    warning('InterneuronPlot:CannotWriteSliceWorkbook', ...
        'Could not write slice audit workbook: %s', err.message);
    writetable(audit.datasetPanel, fullfile(outputFolder, ...
        'interneuron_analyzed_slices_dataset_panel.csv'));
    writetable(audit.datasetMarker, fullfile(outputFolder, ...
        'interneuron_analyzed_slices_marker.csv'));
    writetable(audit.mousePanel, fullfile(outputFolder, ...
        'interneuron_analyzed_slices_mouse_panel.csv'));
    writetable(audit.mouseMarkerRegion, fullfile(outputFolder, ...
        'interneuron_analyzed_slices_mouse_marker_region.csv'));
end
end

function write_excluded_rows_workbook(outputFolder, data)
workbookFile = fullfile(outputFolder, 'interneuron_excluded_rows.xlsx');
excludedSlices = get_import_table(data, 'excludedSlices');
excludedRows = get_import_table(data, 'excludedRows');
try
    if isfile(workbookFile)
        delete(workbookFile);
    end
    writetable(excludedSlices, workbookFile, 'Sheet', 'ExcludedSlices');
    writetable(excludedRows, workbookFile, 'Sheet', 'ExcludedRows');
catch err
    warning('InterneuronPlot:CannotWriteExclusionWorkbook', ...
        'Could not write excluded-row workbook: %s', err.message);
    writetable(excludedSlices, fullfile(outputFolder, ...
        'interneuron_excluded_slices.csv'));
    writetable(excludedRows, fullfile(outputFolder, ...
        'interneuron_excluded_rows.csv'));
end
end

function T = get_import_table(data, fieldName)
if isfield(data, fieldName) && istable(data.(fieldName))
    T = data.(fieldName);
else
    T = table();
end
end

function counts = mouse_slice_counts(T)
if isempty(T) || ~any(strcmp(T.Properties.VariableNames, 'Mouse'))
    counts = nan(0, 1);
    return
end

mice = unique(T.Mouse(isfinite(T.Mouse)));
counts = nan(numel(mice), 1);
for mouseIdx = 1:numel(mice)
    rows = T(T.Mouse == mice(mouseIdx), :);
    if any(strcmp(rows.Properties.VariableNames, 'NAnalyzedSlices'))
        values = rows.NAnalyzedSlices;
    elseif any(strcmp(rows.Properties.VariableNames, 'NInputRows'))
        values = rows.NInputRows;
    else
        values = nan;
    end
    values = unique(values(isfinite(values) & values > 0));
    if ~isempty(values)
        counts(mouseIdx) = max(values);
    end
end
end

function factors = mouse_fov_area_factors(T)
factors = mouse_summary_values(T, 'MeanFovAreaFactor', 1);
end

function factors = mouse_normalization_factors(T)
factors = mouse_summary_values(T, 'CountNormalizationFactor', nan);
if all(~isfinite(factors)) && any(strcmp(T.Properties.VariableNames, 'NAnalyzedSlices'))
    factors = mouse_slice_counts(T);
end
end

function values = mouse_summary_values(T, variableName, defaultValue)
if isempty(T) || ~any(strcmp(T.Properties.VariableNames, 'Mouse'))
    values = nan(0, 1);
    return
end

mice = unique(T.Mouse(isfinite(T.Mouse)));
values = nan(numel(mice), 1);
for mouseIdx = 1:numel(mice)
    rows = T(T.Mouse == mice(mouseIdx), :);
    if any(strcmp(rows.Properties.VariableNames, variableName))
        rowValues = rows.(variableName);
    else
        rowValues = defaultValue;
    end
    rowValues = unique(rowValues(isfinite(rowValues) & rowValues > 0));
    if ~isempty(rowValues)
        values(mouseIdx) = max(rowValues);
    elseif isfinite(defaultValue)
        values(mouseIdx) = defaultValue;
    end
end
end

function label = format_mouse_slice_counts(T)
if isempty(T) || ~any(strcmp(T.Properties.VariableNames, 'Mouse'))
    label = "";
    return
end

mice = unique(T.Mouse(isfinite(T.Mouse)));
parts = strings(numel(mice), 1);
for mouseIdx = 1:numel(mice)
    rows = T(T.Mouse == mice(mouseIdx), :);
    count = mouse_slice_counts(rows);
    if isempty(count)
        count = nan;
    else
        count = count(1);
    end
    parts(mouseIdx) = string(sprintf('%.0f=%s', mice(mouseIdx), ...
        char(format_count_value(count))));
end
label = strjoin(parts, "; ");
end

function label = format_mouse_value_counts(T, variableName)
if isempty(T) || ~any(strcmp(T.Properties.VariableNames, 'Mouse'))
    label = "";
    return
end

mice = unique(T.Mouse(isfinite(T.Mouse)));
parts = strings(numel(mice), 1);
for mouseIdx = 1:numel(mice)
    rows = T(T.Mouse == mice(mouseIdx), :);
    values = mouse_summary_values(rows, variableName, nan);
    if isempty(values)
        value = nan;
    else
        value = values(1);
    end
    parts(mouseIdx) = string(sprintf('%.0f=%s', mice(mouseIdx), ...
        char(format_count_value(value))));
end
label = strjoin(parts, "; ");
end

function label = format_slice_summary(sliceCounts)
sliceCounts = sliceCounts(isfinite(sliceCounts));
if isempty(sliceCounts)
    label = "n/a";
elseif max(sliceCounts) - min(sliceCounts) < 1e-10
    label = format_count_value(sliceCounts(1)) + " each";
else
    label = "mean " + format_count_value(mean(sliceCounts)) + ...
        ", range " + format_count_value(min(sliceCounts)) + "-" + ...
        format_count_value(max(sliceCounts));
end
end

function label = format_factor_summary(factors)
factors = factors(isfinite(factors));
if isempty(factors)
    label = "n/a";
elseif max(factors) - min(factors) < 1e-10
    label = format_count_value(factors(1)) + "x";
else
    label = "mean " + format_count_value(mean(factors)) + ...
        "x, range " + format_count_value(min(factors)) + "-" + ...
        format_count_value(max(factors)) + "x";
end
end

function label = mouse_count_label(nMice)
if nMice == 1
    label = "n = 1 mouse";
else
    label = string(sprintf('n = %.0f mice', nMice));
end
end

function n = count_unique_finite(values)
values = values(isfinite(values));
if isempty(values)
    n = 0;
else
    n = numel(unique(values));
end
end

function countVars = count_variable_names(normalizeCountsBySlices)
if normalizeCountsBySlices
    countVars = struct( ...
        'Active', 'ActiveCountPerSlice', ...
        'TotalRegion', 'TotalCountRegionPerSlice', ...
        'TotalAllRegions', 'TotalCountAllRegionsPerSlice');
else
    countVars = struct( ...
        'Active', 'ActiveCount', ...
        'TotalRegion', 'TotalCountRegion', ...
        'TotalAllRegions', 'TotalCountAllRegions');
end
end

function label = count_unit_label(normalizeCountsBySlices)
if normalizeCountsBySlices
    label = "cells per analyzed slice/stack per mouse (60x-equivalent)";
else
    label = "cells per mouse";
end
end

function formatText = count_value_format(normalizeCountsBySlices)
if normalizeCountsBySlices
    formatText = '%.2f';
else
    formatText = '%.1f';
end
end

function label = format_count_pair(activeValue, totalValue, suffix)
label = format_count_value(activeValue) + "/" + format_count_value(totalValue) + string(suffix);
end

function label = format_count_value(value)
if ~isfinite(value)
    label = "n/a";
elseif abs(value) >= 10
    label = string(sprintf('%.1f', value));
elseif abs(value - round(value)) < 1e-10
    label = string(sprintf('%.0f', value));
else
    label = string(sprintf('%.2f', value));
end
end

function regions = ordered_regions(regionValues)
preferred = ["S Oriens"; "Deep S Pyr"; "Superficial S Pyr"; "S Radiatum"];
regions = strings(0, 1);
for idx = 1:numel(preferred)
    if any(regionValues == preferred(idx))
        regions(end + 1, 1) = preferred(idx); %#ok<AGROW>
    end
end
remaining = unique(regionValues);
remaining = remaining(~ismember(remaining, regions) & remaining ~= "AllRegions");
regions = [regions; remaining(:)];
end

function labels = short_region_labels(regions)
labels = regions;
labels(regions == "S Oriens") = "Oriens";
labels(regions == "Deep S Pyr") = "Deep Pyr";
labels(regions == "Superficial S Pyr") = "Sup Pyr";
labels(regions == "S Radiatum") = "Radiatum";
end

function labels = anatomical_region_labels(regions)
labels = regions;
labels(regions == "S Oriens") = "S. Oriens";
labels(regions == "Deep S Pyr") = "Deep S. Pyr";
labels(regions == "Superficial S Pyr") = "Superficial S. Pyr";
labels(regions == "S Radiatum") = "S. Radiatum";
end

function markers = ordered_markers(markerValues, markerPanel)
panel = string(markerPanel);
if panel == "PV/SOM"
    preferred = ["PV"; "SST/SOM"; "PV+SST/SOM"; "PV-/SST-"];
elseif panel == "SNCG/CCK8"
    preferred = ["SNCG"; "CCK8"; "SNCG+CCK8"; "SNCG-/CCK8-"];
elseif panel == "proCCK/VIP"
    preferred = ["proCCK"; "VIP"; "proCCK+VIP"; "proCCK-/VIP-"];
else
    preferred = unique(markerValues);
end

markers = strings(0, 1);
for idx = 1:numel(preferred)
    if any(markerValues == preferred(idx))
        markers(end + 1, 1) = preferred(idx); %#ok<AGROW>
    end
end
remaining = unique(markerValues);
remaining = remaining(~ismember(remaining, markers));
markers = [markers; remaining(:)];
end

function colors = marker_colors(markers)
fallback = [
    0.20 0.40 0.60
    0.85 0.45 0.12
    0.30 0.55 0.35
    0.45 0.35 0.65
    0.20 0.60 0.60];
colors = nan(numel(markers), 3);
for markerIdx = 1:numel(markers)
    marker = marker_name_from_label(markers(markerIdx));
    if marker == "PV"
        colors(markerIdx, :) = [0.85 0.08 0.08];
    elseif marker == "SST/SOM" || marker == "SOM"
        colors(markerIdx, :) = [0.05 0.22 0.85];
    elseif marker == "PV+SST/SOM" || marker == "SST/SOM+PV" || marker == "PV+SOM"
        colors(markerIdx, :) = [0.50 0.18 0.70];
    elseif marker == "proCCK" || marker == "CCK"
        colors(markerIdx, :) = [0.95 0.72 0.05];
    elseif marker == "VIP"
        colors(markerIdx, :) = [0.05 0.55 0.18];
    elseif marker == "proCCK+VIP" || marker == "VIP+proCCK" || ...
            marker == "VIP+CCK" || marker == "CCK+VIP"
        colors(markerIdx, :) = [0.50 0.28 0.10];
    else
        colors(markerIdx, :) = fallback(mod(markerIdx - 1, size(fallback, 1)) + 1, :);
    end
end
end

function marker = marker_name_from_label(label)
parts = split(string(label), " | ");
marker = strtrim(parts(end));
end

function colorOut = lighten_color(colorIn, amount)
colorOut = colorIn + (1 - colorIn) * amount;
end

function figureInfo = add_figure_info(figureInfo, panel, markerPanel, figureType, fileBase)
newRow = table(panel.Dataset, string(markerPanel), string(figureType), string(fileBase), ...
    'VariableNames', figureInfo.Properties.VariableNames);
figureInfo = [figureInfo; newRow];
end

function save_figure(fig, outputFolder, baseName)
pngFile = fullfile(outputFolder, [baseName '.png']);
pdfFile = fullfile(outputFolder, [baseName '.pdf']);
if exist('exportgraphics', 'file') == 2
    exportgraphics(fig, pngFile, 'Resolution', 300);
    set(fig, 'Renderer', 'painters');
    exportgraphics(fig, pdfFile, 'ContentType', 'vector');
else
    saveas(fig, pngFile);
    set(fig, 'Renderer', 'painters');
    print(fig, pdfFile, '-dpdf', '-painters');
end
end

function write_plot_notes(outputFolder, data, figureInfo, opts, barAxisLimits)
noteFile = fullfile(outputFolder, 'interneuron_plotting_notes.txt');
fid = fopen(noteFile, 'w');
if fid < 0
    warning('InterneuronPlot:CannotWriteNotes', 'Could not write %s.', noteFile);
    return
end
cleanup = onCleanup(@() fclose(fid));
fprintf(fid, 'Interneuron quantification plotting notes\n');
fprintf(fid, 'Workbook: %s\n', char(data.workbookFile));
if isfield(data, 'loadedFrom')
    fprintf(fid, 'Data loaded from: %s\n', char(data.loadedFrom));
end
if isfield(data, 'cacheFile')
    fprintf(fid, 'Cache file: %s\n', char(data.cacheFile));
end
if isfield(data, 'importedAt')
    fprintf(fid, 'Imported at: %s\n', char(data.importedAt));
end
fprintf(fid, 'Counts normalized to 60x-equivalent stacks before mouse averaging: %d\n', ...
    opts.NormalizeCountsBySlices);
fprintf(fid, '\n');
fprintf(fid, 'Datasets are not merged. intrinsic_cFos reports cFos immunostaining; cfos_TRE_mKate reports mKate/RFP tagged cells.\n');
fprintf(fid, 'For cfos_TRE_mKate, marker totals use explicit RFP/mKate-positive plus RFP/mKate-negative blocks where present.\n');
if opts.NormalizeCountsBySlices
    fprintf(fid, 'Count panels show mouse-level cells per analyzed slice/stack, converted to 60x-equivalent counts before averaging across mice.\n');
    fprintf(fid, 'Intrinsic SNCG/CCK8 and proCCK/VIP rows 3:6 and 17:20 use 6 stacks per mouse and a 3x 20x-to-60x field-of-view layer-width normalization, for an 18-fold count normalization.\n');
else
    fprintf(fid, 'Count panels show summed mouse-level cells, averaged across mice.\n');
end
fprintf(fid, 'Answer heatmaps have three vertically stacked panels: mean total cells, mean active/tagged cells, and mean active/tagged fraction, with regions on the y-axis and cell types on the x-axis.\n');
fprintf(fid, 'Selected answer heatmaps include only PV, SOM, PV/SOM, proCCK, and VIP classes and are rendered 20%% wider than the full answer heatmaps.\n');
fprintf(fid, 'Rows/cells marked with an asterisk lack region-level total counts in the workbook.\n');
fprintf(fid, 'For asterisked cells, the total-count panel is blank and the fraction uses active/all-region marker total where available.\n');
fprintf(fid, 'Distribution profile figures ignore active/tagged status and plot total marker cell counts on the same count scale as the bar figures, as mean +/- SEM across mice.\n');
fprintf(fid, 'Bar figures show total bars and active/tagged overlays, faceted by marker panel.\n');
fprintf(fid, 'Bar figures are horizontal anatomical profiles ordered top-to-bottom as S. Oriens, Deep S. Pyr, Superficial S. Pyr, S. Radiatum.\n');
fprintf(fid, 'Bar panel titles include the number of mice for that marker. X-axis labels include the analyzed stacks per mouse and any field-of-view normalization used for count normalization.\n');
fprintf(fid, 'Bar x-axis limits are shared within each dataset/reporter, rounded to compact endpoints, and saved in interneuron_bar_axis_limits.csv.\n');
fprintf(fid, 'In bar labels, an asterisk after the denominator means the denominator is an all-region marker total, not a region-level total.\n\n');
fprintf(fid, 'cfos_TRE_mKate slice groups with no mKate/RFP-positive cells or no positive interneuron-marker cells are excluded before summaries.\n');
fprintf(fid, 'Slice/stack audit workbook: interneuron_analyzed_slices.xlsx\n');
fprintf(fid, 'Excluded-row audit workbook: interneuron_excluded_rows.xlsx\n\n');
fprintf(fid, 'Shared bar x-axis limits:\n');
for rowIdx = 1:height(barAxisLimits)
    fprintf(fid, '- %s | %s: x max %.4g\n', char(barAxisLimits.Dataset(rowIdx)), ...
        char(barAxisLimits.Reporter(rowIdx)), barAxisLimits.BarXAxisMax(rowIdx));
end
fprintf(fid, '\n');
fprintf(fid, 'Generated figure files:\n');
for rowIdx = 1:height(figureInfo)
    fprintf(fid, '- %s | %s | %s: %s\n', char(figureInfo.Dataset(rowIdx)), ...
        char(figureInfo.MarkerPanel(rowIdx)), char(figureInfo.FigureType(rowIdx)), ...
        char(figureInfo.FileBase(rowIdx)));
end
end

function name = sanitize_filename(textValue)
name = regexprep(string(textValue), '[^A-Za-z0-9]+', '_');
name = regexprep(name, '^_+|_+$', '');
end

function folder = default_output_folder()
codeFolder = fileparts(mfilename('fullpath'));
repoRoot = fullfile(codeFolder, '..', '..');
candidateFolders = {
    fullfile(repoRoot, 'Results')
    fullfile(repoRoot, '..', 'Results')};
folder = fullfile(first_existing_folder(candidateFolders), 'Interneuron_quantification_figures');
end

function folder = first_existing_folder(candidateFolders)
folder = candidateFolders{1};
for idx = 1:numel(candidateFolders)
    if exist(candidateFolders{idx}, 'dir')
        folder = candidateFolders{idx};
        return
    end
end
end

function tf = is_name_value_start(value)
if ~(ischar(value) || (isstring(value) && isscalar(value)))
    tf = false;
    return
end
name = lower(char(string(value)));
validNames = {'outputfolder', 'savefigures', 'closefigures', ...
    'includenegativemarkers', 'normalizecountsbyslices', 'cachefile', 'reimportexcel'};
tf = any(strcmp(name, validNames));
end
