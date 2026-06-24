function out = plot_CB_quantification_figures(workbookFile, varargin)
%PLOT_CB_QUANTIFICATION_FIGURES First-pass CB/RFP quantification figures.

if nargin < 1
    workbookFile = [];
elseif is_name_value_start(workbookFile)
    varargin = [{workbookFile}, varargin];
    workbookFile = [];
end

parser = inputParser;
addParameter(parser, 'OutputFolder', default_output_folder(), ...
    @(x) ischar(x) || isstring(x));
addParameter(parser, 'SideFilter', 'all', @(x) ischar(x) || isstring(x));
addParameter(parser, 'SaveFigures', true, @(x) islogical(x) && isscalar(x));
addParameter(parser, 'CloseFigures', false, @(x) islogical(x) && isscalar(x));
addParameter(parser, 'CacheFile', '', @(x) ischar(x) || isstring(x));
addParameter(parser, 'ReimportExcel', false, @(x) islogical(x) && isscalar(x));
parse(parser, varargin{:});
opts = parser.Results;
opts.OutputFolder = char(opts.OutputFolder);
opts.SideFilter = string(opts.SideFilter);

data = import_CB_quantifications(workbookFile, 'IncludeZeroCounts', true, ...
    'CacheFile', opts.CacheFile, 'ReimportExcel', opts.ReimportExcel);
longTable = data.long;
longTable = longTable(ismember(longTable.Position, ["Superficial", "Deep"]), :);

if lower(opts.SideFilter) ~= "all"
    longTable = longTable(strcmpi(longTable.Side, opts.SideFilter), :);
end
if isempty(longTable)
    error('CBPlot:NoData', 'No CB quantification rows matched SideFilter=%s.', char(opts.SideFilter));
end

summaryTable = summarize_layer_status(longTable);
pairedTable = make_paired_layer_table(summaryTable);
cbLayerAllocationTable = make_layer_allocation_table(summaryTable);
[layerAllocationTable, layerAllocationSource] = choose_layer_allocation_table(data, cbLayerAllocationTable, opts);
statusOverallTable = summarize_status_overall(summaryTable);
statusCompositeTable = summarize_status_composite_overall(statusOverallTable);
layerCompositeWithinTable = make_layer_composite_within_table(summaryTable);
layerCompositeGlobalTable = make_layer_composite_global_table(summaryTable);

if opts.SaveFigures && ~exist(opts.OutputFolder, 'dir')
    mkdir(opts.OutputFolder);
end

graphValuesWorkbook = "";
prismWorkbook = "";

figHandles = gobjects(8, 1);
figHandles(1) = plot_layer_allocation_paired(layerAllocationTable, opts);
figHandles(2) = plot_status_assignment_overall(statusOverallTable, opts);
figHandles(3) = plot_status_assignment_composite_column(statusCompositeTable, opts);
figHandles(4) = plot_cb_positive_paired(pairedTable, opts);
figHandles(5) = plot_status_split_paired(pairedTable, opts);
figHandles(6) = plot_stacked_composition(summaryTable, opts);
figHandles(7) = plot_layer_status_heatmap(layerCompositeWithinTable, opts);
figHandles(8) = plot_layer_status_global_heatmap(layerCompositeGlobalTable, opts);

if opts.SaveFigures
    writetable(summaryTable, fullfile(opts.OutputFolder, 'cb_layer_status_summary.csv'));
    writetable(pairedTable, fullfile(opts.OutputFolder, 'cb_layer_paired_percentages.csv'));
    writetable(layerAllocationTable, fullfile(opts.OutputFolder, 'cb_mkate_layer_allocation_summary.csv'));
    writetable(statusOverallTable, fullfile(opts.OutputFolder, 'cb_status_overall_assignment_summary.csv'));
    writetable(statusCompositeTable, fullfile(opts.OutputFolder, 'cb_status_composite_assignment_summary.csv'));
    writetable(layerCompositeWithinTable, fullfile(opts.OutputFolder, 'cb_layer_composite_within_layer_heatmap_summary.csv'));
    writetable(layerCompositeGlobalTable, fullfile(opts.OutputFolder, 'cb_layer_composite_global_heatmap_summary.csv'));
    save_figure(figHandles(1), opts.OutputFolder, 'cb_mkate_layer_allocation_paired');
    save_figure(figHandles(2), opts.OutputFolder, 'cb_status_overall_assignment');
    save_figure(figHandles(3), opts.OutputFolder, 'cb_status_composite_assignment_column');
    save_figure(figHandles(4), opts.OutputFolder, 'cb_positive_deep_vs_superficial_paired');
    save_figure(figHandles(5), opts.OutputFolder, 'cb_status_deep_vs_superficial_paired');
    save_figure(figHandles(6), opts.OutputFolder, 'cb_status_stacked_composition_by_mouse_layer');
    save_figure(figHandles(7), opts.OutputFolder, 'cb_status_layer_percent_heatmap');
    save_figure(figHandles(8), opts.OutputFolder, 'cb_status_layer_global_percent_heatmap');
    [graphValuesWorkbook, prismWorkbook] = export_graph_value_workbooks(opts.OutputFolder, ...
        data, opts, layerAllocationSource, layerAllocationTable, statusOverallTable, statusCompositeTable, ...
        pairedTable, summaryTable, layerCompositeWithinTable, layerCompositeGlobalTable);
    write_panel_guide(opts.OutputFolder, data, opts, layerAllocationSource, height(layerAllocationTable));
    write_overview_html(opts.OutputFolder, data, opts, layerAllocationSource, height(layerAllocationTable));
end

if opts.CloseFigures
    close(figHandles(ishandle(figHandles)));
end

out = struct();
out.data = data;
out.layerStatusSummary = summaryTable;
out.pairedLayerPercentages = pairedTable;
out.layerAllocation = layerAllocationTable;
out.cbLayerAllocation = cbLayerAllocationTable;
out.layerAllocationSource = layerAllocationSource;
out.statusOverallAssignment = statusOverallTable;
out.statusCompositeAssignment = statusCompositeTable;
out.layerCompositeWithin = layerCompositeWithinTable;
out.layerCompositeGlobal = layerCompositeGlobalTable;
out.figures = figHandles;
out.outputFolder = string(opts.OutputFolder);
out.overviewFile = string(fullfile(opts.OutputFolder, 'overview.html'));
out.graphValuesWorkbook = string(graphValuesWorkbook);
out.prismWorkbook = string(prismWorkbook);
out.sideFilter = opts.SideFilter;
end

function summaryTable = summarize_layer_status(longTable)
mice = unique(longTable.Mouse);
mice = mice(~isnan(mice));
layers = ["Superficial"; "Deep"];
statuses = ["CB+"; "CB-"; "CB+/-"];

nRows = numel(mice) * numel(layers) * numel(statuses);
Mouse = nan(nRows, 1);
Layer = strings(nRows, 1);
CbStatus = strings(nRows, 1);
Count = nan(nRows, 1);
TotalLayerCount = nan(nRows, 1);
PercentOfLayer = nan(nRows, 1);

k = 0;
for mouseIdx = 1:numel(mice)
    for layerIdx = 1:numel(layers)
        layerMask = longTable.Mouse == mice(mouseIdx) & longTable.Position == layers(layerIdx);
        layerTotal = sum(longTable.Count(layerMask));
        for statusIdx = 1:numel(statuses)
            k = k + 1;
            statusMask = layerMask & longTable.CbStatus == statuses(statusIdx);
            Mouse(k) = mice(mouseIdx);
            Layer(k) = layers(layerIdx);
            CbStatus(k) = statuses(statusIdx);
            Count(k) = sum(longTable.Count(statusMask));
            TotalLayerCount(k) = layerTotal;
            PercentOfLayer(k) = safe_percent(Count(k), layerTotal);
        end
    end
end

summaryTable = table(Mouse, Layer, CbStatus, Count, TotalLayerCount, PercentOfLayer);
end

function pairedTable = make_paired_layer_table(summaryTable)
mice = unique(summaryTable.Mouse);
statuses = ["CB+"; "CB-"; "CB+/-"];
nRows = numel(mice) * numel(statuses);

Mouse = nan(nRows, 1);
CbStatus = strings(nRows, 1);
SuperficialPercent = nan(nRows, 1);
DeepPercent = nan(nRows, 1);
DeepMinusSuperficialPercent = nan(nRows, 1);
SuperficialCount = nan(nRows, 1);
DeepCount = nan(nRows, 1);
SuperficialTotal = nan(nRows, 1);
DeepTotal = nan(nRows, 1);

k = 0;
for mouseIdx = 1:numel(mice)
    for statusIdx = 1:numel(statuses)
        k = k + 1;
        Mouse(k) = mice(mouseIdx);
        CbStatus(k) = statuses(statusIdx);
        supRow = summaryTable.Mouse == mice(mouseIdx) & ...
            summaryTable.Layer == "Superficial" & summaryTable.CbStatus == statuses(statusIdx);
        deepRow = summaryTable.Mouse == mice(mouseIdx) & ...
            summaryTable.Layer == "Deep" & summaryTable.CbStatus == statuses(statusIdx);
        SuperficialPercent(k) = first_value(summaryTable.PercentOfLayer(supRow));
        DeepPercent(k) = first_value(summaryTable.PercentOfLayer(deepRow));
        DeepMinusSuperficialPercent(k) = DeepPercent(k) - SuperficialPercent(k);
        SuperficialCount(k) = first_value(summaryTable.Count(supRow));
        DeepCount(k) = first_value(summaryTable.Count(deepRow));
        SuperficialTotal(k) = first_value(summaryTable.TotalLayerCount(supRow));
        DeepTotal(k) = first_value(summaryTable.TotalLayerCount(deepRow));
    end
end

pairedTable = table(Mouse, CbStatus, SuperficialPercent, DeepPercent, ...
    DeepMinusSuperficialPercent, SuperficialCount, DeepCount, SuperficialTotal, DeepTotal);
end

function layerAllocationTable = make_layer_allocation_table(summaryTable)
mice = unique(summaryTable.Mouse);

Mouse = nan(numel(mice), 1);
SuperficialCount = nan(numel(mice), 1);
DeepCount = nan(numel(mice), 1);
InLayerTotal = nan(numel(mice), 1);
SuperficialPercentOfInLayer = nan(numel(mice), 1);
DeepPercentOfInLayer = nan(numel(mice), 1);
DeepMinusSuperficialPercent = nan(numel(mice), 1);

for mouseIdx = 1:numel(mice)
    Mouse(mouseIdx) = mice(mouseIdx);
    supRow = summaryTable.Mouse == mice(mouseIdx) & summaryTable.Layer == "Superficial";
    deepRow = summaryTable.Mouse == mice(mouseIdx) & summaryTable.Layer == "Deep";
    SuperficialCount(mouseIdx) = first_value(summaryTable.TotalLayerCount(supRow));
    DeepCount(mouseIdx) = first_value(summaryTable.TotalLayerCount(deepRow));
    InLayerTotal(mouseIdx) = SuperficialCount(mouseIdx) + DeepCount(mouseIdx);
    SuperficialPercentOfInLayer(mouseIdx) = safe_percent(SuperficialCount(mouseIdx), InLayerTotal(mouseIdx));
    DeepPercentOfInLayer(mouseIdx) = safe_percent(DeepCount(mouseIdx), InLayerTotal(mouseIdx));
    DeepMinusSuperficialPercent(mouseIdx) = DeepPercentOfInLayer(mouseIdx) - SuperficialPercentOfInLayer(mouseIdx);
end

layerAllocationTable = table(Mouse, SuperficialCount, DeepCount, InLayerTotal, ...
    SuperficialPercentOfInLayer, DeepPercentOfInLayer, DeepMinusSuperficialPercent);
end

function [layerAllocationTable, layerAllocationSource] = choose_layer_allocation_table(data, cbLayerAllocationTable, opts)
layerAllocationTable = cbLayerAllocationTable;
layerAllocationSource = sprintf('%s sheet', char(data.sheetName));

if lower(opts.SideFilter) ~= "all"
    return
end

if isfield(data, 'mKate') && isfield(data.mKate, 'layerAllocation') && ...
        ~isempty(data.mKate.layerAllocation) && height(data.mKate.layerAllocation) > 0
    layerAllocationTable = combine_layer_allocation_tables(cbLayerAllocationTable, data.mKate.layerAllocation);
    layerAllocationSource = sprintf('%s + %s sheets, matched by Mouse', ...
        char(data.sheetName), char(data.mKate.sheetName));
end
end

function layerAllocationTable = combine_layer_allocation_tables(cbLayerAllocationTable, mKateLayerAllocationTable)
mice = unique([cbLayerAllocationTable.Mouse; mKateLayerAllocationTable.Mouse]);
mice = mice(~isnan(mice));

Mouse = mice(:);
CbSuperficialCount = nan(numel(mice), 1);
CbDeepCount = nan(numel(mice), 1);
MKateSuperficialCount = nan(numel(mice), 1);
MKateDeepCount = nan(numel(mice), 1);
SuperficialCount = nan(numel(mice), 1);
DeepCount = nan(numel(mice), 1);
InLayerTotal = nan(numel(mice), 1);
SuperficialPercentOfInLayer = nan(numel(mice), 1);
DeepPercentOfInLayer = nan(numel(mice), 1);
DeepMinusSuperficialPercent = nan(numel(mice), 1);
Source = strings(numel(mice), 1);

for mouseIdx = 1:numel(mice)
    mouse = mice(mouseIdx);
    cbRows = cbLayerAllocationTable.Mouse == mouse;
    mKateRows = mKateLayerAllocationTable.Mouse == mouse;

    cbSuperficial = sum_table_values(cbLayerAllocationTable, cbRows, 'SuperficialCount');
    cbDeep = sum_table_values(cbLayerAllocationTable, cbRows, 'DeepCount');
    mKateSuperficial = sum_table_values(mKateLayerAllocationTable, mKateRows, 'SuperficialCount');
    mKateDeep = sum_table_values(mKateLayerAllocationTable, mKateRows, 'DeepCount');

    CbSuperficialCount(mouseIdx) = cbSuperficial;
    CbDeepCount(mouseIdx) = cbDeep;
    MKateSuperficialCount(mouseIdx) = mKateSuperficial;
    MKateDeepCount(mouseIdx) = mKateDeep;
    SuperficialCount(mouseIdx) = cbSuperficial + mKateSuperficial;
    DeepCount(mouseIdx) = cbDeep + mKateDeep;
    InLayerTotal(mouseIdx) = SuperficialCount(mouseIdx) + DeepCount(mouseIdx);
    SuperficialPercentOfInLayer(mouseIdx) = safe_percent(SuperficialCount(mouseIdx), InLayerTotal(mouseIdx));
    DeepPercentOfInLayer(mouseIdx) = safe_percent(DeepCount(mouseIdx), InLayerTotal(mouseIdx));
    DeepMinusSuperficialPercent(mouseIdx) = DeepPercentOfInLayer(mouseIdx) - SuperficialPercentOfInLayer(mouseIdx);

    if any(cbRows) && any(mKateRows)
        Source(mouseIdx) = "CB+RFP + MKate";
    elseif any(cbRows)
        Source(mouseIdx) = "CB+RFP";
    else
        Source(mouseIdx) = "MKate";
    end
end

layerAllocationTable = table(Mouse, SuperficialCount, DeepCount, InLayerTotal, ...
    SuperficialPercentOfInLayer, DeepPercentOfInLayer, DeepMinusSuperficialPercent, ...
    CbSuperficialCount, CbDeepCount, MKateSuperficialCount, MKateDeepCount, Source);
end

function value = sum_table_values(sourceTable, rowMask, variableName)
if ~any(rowMask)
    value = 0;
else
    values = sourceTable.(variableName)(rowMask);
    values = values(isfinite(values));
    value = sum(values);
end
end

function statusOverallTable = summarize_status_overall(summaryTable)
mice = unique(summaryTable.Mouse);
statuses = ["CB+"; "CB-"; "CB+/-"];
nRows = numel(mice) * numel(statuses);

Mouse = nan(nRows, 1);
CbStatus = strings(nRows, 1);
Count = nan(nRows, 1);
InLayerTotal = nan(nRows, 1);
PercentOfInLayerMKate = nan(nRows, 1);

k = 0;
for mouseIdx = 1:numel(mice)
    mouseRows = summaryTable.Mouse == mice(mouseIdx);
    mouseTotal = sum(summaryTable.Count(mouseRows));
    for statusIdx = 1:numel(statuses)
        k = k + 1;
        statusRows = mouseRows & summaryTable.CbStatus == statuses(statusIdx);
        Mouse(k) = mice(mouseIdx);
        CbStatus(k) = statuses(statusIdx);
        Count(k) = sum(summaryTable.Count(statusRows));
        InLayerTotal(k) = mouseTotal;
        PercentOfInLayerMKate(k) = safe_percent(Count(k), mouseTotal);
    end
end

statusOverallTable = table(Mouse, CbStatus, Count, InLayerTotal, PercentOfInLayerMKate);
end

function statusCompositeTable = summarize_status_composite_overall(statusOverallTable)
mice = unique(statusOverallTable.Mouse);

Mouse = nan(numel(mice), 1);
CbNegativeCount = nan(numel(mice), 1);
CbPositiveCount = nan(numel(mice), 1);
CbIntermediateCount = nan(numel(mice), 1);
CbPositiveIntermediateCount = nan(numel(mice), 1);
InLayerTotal = nan(numel(mice), 1);
CbNegativePercent = nan(numel(mice), 1);
CbPositivePercent = nan(numel(mice), 1);
CbIntermediatePercent = nan(numel(mice), 1);
CbPositiveIntermediatePercent = nan(numel(mice), 1);

for mouseIdx = 1:numel(mice)
    mouse = mice(mouseIdx);
    Mouse(mouseIdx) = mouse;
    cbNegativeRow = statusOverallTable.Mouse == mouse & statusOverallTable.CbStatus == "CB-";
    cbPositiveRow = statusOverallTable.Mouse == mouse & statusOverallTable.CbStatus == "CB+";
    cbIntermediateRow = statusOverallTable.Mouse == mouse & statusOverallTable.CbStatus == "CB+/-";

    CbNegativeCount(mouseIdx) = first_value(statusOverallTable.Count(cbNegativeRow));
    CbPositiveCount(mouseIdx) = first_value(statusOverallTable.Count(cbPositiveRow));
    CbIntermediateCount(mouseIdx) = first_value(statusOverallTable.Count(cbIntermediateRow));
    CbPositiveIntermediateCount(mouseIdx) = CbPositiveCount(mouseIdx) + CbIntermediateCount(mouseIdx);
    InLayerTotal(mouseIdx) = first_value(statusOverallTable.InLayerTotal(cbNegativeRow));
    CbNegativePercent(mouseIdx) = safe_percent(CbNegativeCount(mouseIdx), InLayerTotal(mouseIdx));
    CbPositivePercent(mouseIdx) = safe_percent(CbPositiveCount(mouseIdx), InLayerTotal(mouseIdx));
    CbIntermediatePercent(mouseIdx) = safe_percent(CbIntermediateCount(mouseIdx), InLayerTotal(mouseIdx));
    CbPositiveIntermediatePercent(mouseIdx) = safe_percent(CbPositiveIntermediateCount(mouseIdx), InLayerTotal(mouseIdx));
end

statusCompositeTable = table(Mouse, CbNegativeCount, CbPositiveCount, CbIntermediateCount, ...
    CbPositiveIntermediateCount, InLayerTotal, CbNegativePercent, CbPositivePercent, ...
    CbIntermediatePercent, CbPositiveIntermediatePercent);
end

function layerCompositeWithinTable = make_layer_composite_within_table(summaryTable)
layerCompositeWithinTable = make_layer_composite_table(summaryTable, "within_layer");
end

function layerCompositeGlobalTable = make_layer_composite_global_table(summaryTable)
layerCompositeGlobalTable = make_layer_composite_table(summaryTable, "superficial_deep_total");
end

function compositeTable = make_layer_composite_table(summaryTable, denominatorMode)
mice = unique(summaryTable.Mouse);
layers = ["Superficial"; "Deep"];
groups = ["CB+/intermediate"; "CB-"];
nRows = numel(mice) * numel(layers) * numel(groups);

Mouse = nan(nRows, 1);
Layer = strings(nRows, 1);
CbGroup = strings(nRows, 1);
Count = nan(nRows, 1);
Denominator = nan(nRows, 1);
Percent = nan(nRows, 1);

k = 0;
for mouseIdx = 1:numel(mice)
    mouse = mice(mouseIdx);
    mouseRows = summaryTable.Mouse == mouse;
    globalDenominator = sum(summaryTable.Count(mouseRows));
    for layerIdx = 1:numel(layers)
        layer = layers(layerIdx);
        layerRows = mouseRows & summaryTable.Layer == layer;
        layerDenominator = first_value(summaryTable.TotalLayerCount(layerRows));
        if denominatorMode == "within_layer"
            denominator = layerDenominator;
        else
            denominator = globalDenominator;
        end

        cbPositiveCount = get_layer_status_count(summaryTable, mouse, layer, "CB+");
        cbIntermediateCount = get_layer_status_count(summaryTable, mouse, layer, "CB+/-");
        cbNegativeCount = get_layer_status_count(summaryTable, mouse, layer, "CB-");
        counts = [cbPositiveCount + cbIntermediateCount, cbNegativeCount];

        for groupIdx = 1:numel(groups)
            k = k + 1;
            Mouse(k) = mouse;
            Layer(k) = layer;
            CbGroup(k) = groups(groupIdx);
            Count(k) = counts(groupIdx);
            Denominator(k) = denominator;
            Percent(k) = safe_percent(Count(k), Denominator(k));
        end
    end
end

compositeTable = table(Mouse, Layer, CbGroup, Count, Denominator, Percent);
end

function count = get_layer_status_count(summaryTable, mouse, layer, cbStatus)
row = summaryTable.Mouse == mouse & summaryTable.Layer == layer & summaryTable.CbStatus == cbStatus;
count = first_value(summaryTable.Count(row));
end

function fig = plot_layer_allocation_paired(layerAllocationTable, opts)
fig = figure('Color', 'w', 'Units', 'centimeters', 'Position', [4 4 8 8]);
ax = axes(fig);
hold(ax, 'on');
layerColors = [0.05 0.45 0.55; 0.78 0.36 0.12];

for rowIdx = 1:height(layerAllocationTable)
    y = [layerAllocationTable.SuperficialPercentOfInLayer(rowIdx), ...
        layerAllocationTable.DeepPercentOfInLayer(rowIdx)];
    plot(ax, [1 2], y, '-', 'Color', [0.72 0.72 0.72], 'LineWidth', 0.8);
    scatter(ax, 1, y(1), 36, layerColors(1, :), 'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 0.3);
    scatter(ax, 2, y(2), 36, layerColors(2, :), 'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 0.3);
end

plot_mean_sem(ax, 1, layerAllocationTable.SuperficialPercentOfInLayer);
plot_mean_sem(ax, 2, layerAllocationTable.DeepPercentOfInLayer);
set(ax, 'XLim', [0.65 2.35], 'XTick', [1 2], 'XTickLabel', {'Superficial', 'Deep'}, ...
    'TickDir', 'out');
ylabel(ax, '% of in-layer mKate/RFP cells');
ylim(ax, [0 100]);
title(ax, {title_with_filter('Engram allocation across CA1 layers', opts), ...
    '% of mKate/RFP cells in superficial vs deep CA1'}, 'Interpreter', 'none');
box(ax, 'off');
hold(ax, 'off');
end

function fig = plot_status_assignment_overall(statusOverallTable, opts)
fig = figure('Color', 'w', 'Units', 'centimeters', 'Position', [4 4 10 8]);
ax = axes(fig);
hold(ax, 'on');
statuses = ["CB+"; "CB-"; "CB+/-"];
colors = cb_status_colors();
mice = unique(statusOverallTable.Mouse);

for mouseIdx = 1:numel(mice)
    mouse = mice(mouseIdx);
    y = nan(1, numel(statuses));
    for statusIdx = 1:numel(statuses)
        row = statusOverallTable.Mouse == mouse & statusOverallTable.CbStatus == statuses(statusIdx);
        y(statusIdx) = first_value(statusOverallTable.PercentOfInLayerMKate(row));
    end
    plot(ax, 1:numel(statuses), y, '-', 'Color', [0.75 0.75 0.75], 'LineWidth', 0.8);
    for statusIdx = 1:numel(statuses)
        scatter(ax, statusIdx, y(statusIdx), 36, colors(statusIdx, :), ...
            'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 0.3);
    end
end

for statusIdx = 1:numel(statuses)
    values = statusOverallTable.PercentOfInLayerMKate(statusOverallTable.CbStatus == statuses(statusIdx));
    plot_mean_sem(ax, statusIdx, values);
end

set(ax, 'XLim', [0.65 3.35], 'XTick', 1:3, 'XTickLabel', cellstr(statuses), 'TickDir', 'out');
ylabel(ax, '% of in-layer mKate/RFP cells');
ylim(ax, [0 100]);
title(ax, {title_with_filter('Engram assignment to CB status', opts), ...
    'Overall superficial+deep CA1 mKate/RFP cells'}, 'Interpreter', 'none');
box(ax, 'off');
hold(ax, 'off');
end

function fig = plot_status_assignment_composite_column(statusCompositeTable, opts)
fig = figure('Color', 'w', 'Units', 'centimeters', 'Position', [4 4 10 8]);
ax = axes(fig);
hold(ax, 'on');
colors = cb_status_colors();

meanCbNegative = mean_omit_nan(statusCompositeTable.CbNegativePercent);
meanCbPositive = mean_omit_nan(statusCompositeTable.CbPositivePercent);
meanCbIntermediate = mean_omit_nan(statusCompositeTable.CbIntermediatePercent);
barValues = [
    meanCbNegative, 0, 0
    0, meanCbPositive, meanCbIntermediate];

b = bar(ax, barValues, 'stacked', 'LineWidth', 0.5);
b(1).FaceColor = colors(2, :);
b(2).FaceColor = colors(1, :);
b(3).FaceColor = colors(3, :);

for rowIdx = 1:height(statusCompositeTable)
    y = [statusCompositeTable.CbNegativePercent(rowIdx), ...
        statusCompositeTable.CbPositiveIntermediatePercent(rowIdx)];
    plot(ax, [1 2], y, '-', 'Color', [0.75 0.75 0.75], 'LineWidth', 0.8, ...
        'HandleVisibility', 'off');
    scatter(ax, 1, y(1), 36, colors(2, :), 'filled', ...
        'MarkerEdgeColor', 'k', 'LineWidth', 0.3, 'HandleVisibility', 'off');
    scatter(ax, 2, y(2), 36, [0.15 0.15 0.15], 'filled', ...
        'MarkerEdgeColor', 'k', 'LineWidth', 0.3, 'HandleVisibility', 'off');
end

set(ax, 'XLim', [0.5 2.5], 'XTick', [1 2], ...
    'XTickLabel', {'CB-', 'CB+ + intermediate'}, 'TickDir', 'out');
ylabel(ax, '% of in-layer mKate/RFP cells');
ylim(ax, [0 100]);
legend(ax, b, {'CB-', 'CB+', 'CB+/-'}, 'Location', 'eastoutside');
title(ax, {title_with_filter('Composite CB assignment', opts), ...
    'Column 2 stacks CB+ and CB+/- intermediate cells'}, 'Interpreter', 'none');
box(ax, 'off');
hold(ax, 'off');
end

function fig = plot_cb_positive_paired(pairedTable, opts)
fig = figure('Color', 'w', 'Units', 'centimeters', 'Position', [4 4 8 8]);
ax = axes(fig);
plot_paired_layer_values(ax, pairedTable, "CB+", '% CB+ among mKate/RFP cells', opts);
title(ax, {title_with_filter('CB+ fraction by CA1 layer', opts), ...
    'Per layer: % of mKate/RFP cells that are CB+'}, 'Interpreter', 'none');
end

function fig = plot_status_split_paired(pairedTable, opts)
fig = figure('Color', 'w', 'Units', 'centimeters', 'Position', [4 4 18 7]);
tiledlayout(fig, 1, 3, 'Padding', 'compact', 'TileSpacing', 'compact');
statuses = ["CB+"; "CB-"; "CB+/-"];
for statusIdx = 1:numel(statuses)
    ax = nexttile;
    plot_paired_layer_values(ax, pairedTable, statuses(statusIdx), ...
        sprintf('%% %s among mKate/RFP cells', statuses(statusIdx)), opts);
    title(ax, statuses(statusIdx), 'Interpreter', 'none');
end
sgtitle(title_with_filter('Per-layer mKate/RFP assignment to CB status', opts), 'Interpreter', 'none');
end

function fig = plot_stacked_composition(summaryTable, opts)
statuses = ["CB+"; "CB-"; "CB+/-"];
layers = ["Superficial"; "Deep"];
mice = unique(summaryTable.Mouse);
barValues = nan(numel(mice) * numel(layers), numel(statuses));
barTotals = nan(numel(mice) * numel(layers), 1);
tickLabels = cell(numel(mice) * numel(layers), 1);

k = 0;
for mouseIdx = 1:numel(mice)
    for layerIdx = 1:numel(layers)
        k = k + 1;
        for statusIdx = 1:numel(statuses)
            row = summaryTable.Mouse == mice(mouseIdx) & ...
                summaryTable.Layer == layers(layerIdx) & summaryTable.CbStatus == statuses(statusIdx);
            barValues(k, statusIdx) = first_value(summaryTable.PercentOfLayer(row));
            barTotals(k) = first_value(summaryTable.TotalLayerCount(row));
        end
        tickLabels{k} = sprintf('%g\n%s', mice(mouseIdx), layer_short_name(layers(layerIdx)));
    end
end

fig = figure('Color', 'w', 'Units', 'centimeters', 'Position', [4 4 17 8]);
ax = axes(fig);
b = bar(ax, barValues, 'stacked', 'LineWidth', 0.4);
colors = cb_status_colors();
for statusIdx = 1:numel(statuses)
    b(statusIdx).FaceColor = colors(statusIdx, :);
end
set(ax, 'XTick', 1:numel(tickLabels), 'XTickLabel', tickLabels, 'TickDir', 'out');
ylabel(ax, '% of mKate/RFP cells in layer');
ylim(ax, [0 108]);
title(ax, {title_with_filter('CB status composition by mouse and layer', opts), ...
    'Each bar sums to 100% of mKate/RFP cells in that mouse/layer'}, 'Interpreter', 'none');
legend(ax, b, cellstr(statuses), 'Location', 'eastoutside');
box(ax, 'off');
for idx = 1:numel(barTotals)
    text(ax, idx, 102, sprintf('n=%g', barTotals(idx)), ...
        'HorizontalAlignment', 'center', 'FontSize', 7);
end
for mouseIdx = 1:(numel(mice) - 1)
    line(ax, [mouseIdx * 2 + 0.5, mouseIdx * 2 + 0.5], [0 100], 'LineStyle', ':', ...
        'Color', [0.75 0.75 0.75], 'LineWidth', 0.75, 'HandleVisibility', 'off');
end
end

function fig = plot_layer_status_heatmap(layerCompositeWithinTable, opts)
[heatValues, mice] = make_composite_heatmap_matrix(layerCompositeWithinTable);
fig = plot_composite_heatmap(heatValues, mice, ...
    {'S CB+/int', 'S CB-', 'D CB+/int', 'D CB-'}, ...
    '% of mKate/RFP cells within mouse/layer', ...
    'S = superficial, D = deep; each layer pair sums to 100% per mouse', ...
    {title_with_filter('CB status by layer heatmap', opts), ...
    'CB+ and CB+/- intermediate combined; each S or D pair sums to 100%'});
end

function fig = plot_layer_status_global_heatmap(layerCompositeGlobalTable, opts)
[heatValues, mice] = make_composite_heatmap_matrix(layerCompositeGlobalTable);
fig = plot_composite_heatmap(heatValues, mice, ...
    {'S CB+/int', 'S CB-', 'D CB+/int', 'D CB-'}, ...
    '% of superficial+deep mKate/RFP cells', ...
    'All four columns sum to 100% per mouse', ...
    {title_with_filter('Global layer/status heatmap', opts), ...
    'CB+ and CB+/- intermediate combined; all four columns sum to 100%'});
end

function [heatValues, mice] = make_composite_heatmap_matrix(compositeTable)
mice = unique(compositeTable.Mouse);
heatValues = nan(numel(mice), 4);
layers = ["Superficial"; "Deep"];
groups = ["CB+/intermediate"; "CB-"];

for mouseIdx = 1:numel(mice)
    colIdx = 0;
    for layerIdx = 1:numel(layers)
        for groupIdx = 1:numel(groups)
            colIdx = colIdx + 1;
            row = compositeTable.Mouse == mice(mouseIdx) & ...
                compositeTable.Layer == layers(layerIdx) & compositeTable.CbGroup == groups(groupIdx);
            heatValues(mouseIdx, colIdx) = first_value(compositeTable.Percent(row));
        end
    end
end
end

function fig = plot_composite_heatmap(heatValues, mice, xTickLabels, colorbarLabel, xLabelText, titleText)
fig = figure('Color', 'w', 'Units', 'centimeters', 'Position', [4 4 12 7]);
ax = axes(fig);
imagesc(ax, heatValues);
colormap(ax, parula);
axes(ax);
caxis([0 100]);
colorHandle = colorbar;
ylabel(colorHandle, colorbarLabel);
set(ax, 'XTick', 1:numel(xTickLabels), 'XTickLabel', xTickLabels, ...
    'YTick', 1:numel(mice), 'YTickLabel', cellstr(string(mice)), 'TickDir', 'out');
xtickangle(ax, 45);
xlabel(ax, xLabelText);
ylabel(ax, 'Mouse');
title(ax, titleText, 'Interpreter', 'none');

for rowIdx = 1:size(heatValues, 1)
    for colIdx = 1:size(heatValues, 2)
        value = heatValues(rowIdx, colIdx);
        if isnan(value)
            label = 'NA';
        else
            label = sprintf('%.1f', value);
        end
        textColor = [1 1 1];
        if isnan(value) || value < 55
            textColor = [0.1 0.1 0.1];
        end
        text(ax, colIdx, rowIdx, label, 'HorizontalAlignment', 'center', ...
            'FontSize', 8, 'Color', textColor);
    end
end
end

function plot_paired_layer_values(ax, pairedTable, statusName, yLabelText, opts)
rowMask = pairedTable.CbStatus == statusName;
T = pairedTable(rowMask, :);
hold(ax, 'on');
layerColors = [0.05 0.45 0.55; 0.78 0.36 0.12];

for rowIdx = 1:height(T)
    y = [T.SuperficialPercent(rowIdx), T.DeepPercent(rowIdx)];
    if all(isfinite(y))
        plot(ax, [1 2], y, '-', 'Color', [0.72 0.72 0.72], 'LineWidth', 0.8);
    end
    if isfinite(y(1))
        scatter(ax, 1, y(1), 36, layerColors(1, :), 'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 0.3);
    end
    if isfinite(y(2))
        scatter(ax, 2, y(2), 36, layerColors(2, :), 'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 0.3);
    end
end

plot_mean_sem(ax, 1, T.SuperficialPercent);
plot_mean_sem(ax, 2, T.DeepPercent);
set(ax, 'XLim', [0.65 2.35], 'XTick', [1 2], 'XTickLabel', {'Superficial', 'Deep'}, ...
    'TickDir', 'out');
ylabel(ax, yLabelText, 'Interpreter', 'none');
finiteValues = [T.SuperficialPercent; T.DeepPercent];
finiteValues = finiteValues(isfinite(finiteValues));
if isempty(finiteValues)
    yMax = 100;
else
    yMax = max(100, max(finiteValues) * 1.12);
end
ylim(ax, [0 yMax]);
box(ax, 'off');
hold(ax, 'off');
end

function plot_mean_sem(ax, x, values)
values = values(isfinite(values));
if isempty(values)
    return
end
meanValue = mean(values);
semValue = std(values) / sqrt(numel(values));
errorbar(ax, x, meanValue, semValue, 'ko', 'MarkerFaceColor', 'w', ...
    'MarkerSize', 6, 'LineWidth', 1.2, 'CapSize', 8);
end

function value = mean_omit_nan(values)
values = values(isfinite(values));
if isempty(values)
    value = nan;
else
    value = mean(values);
end
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

function [graphValuesWorkbook, prismWorkbook] = export_graph_value_workbooks(outputFolder, ...
    data, opts, layerAllocationSource, layerAllocationTable, statusOverallTable, statusCompositeTable, ...
    pairedTable, summaryTable, layerCompositeWithinTable, layerCompositeGlobalTable)

graphValuesWorkbook = fullfile(outputFolder, 'cb_graph_values_for_verification.xlsx');
prismWorkbook = fullfile(outputFolder, 'cb_prism_ready_tables.xlsx');

reset_output_file(graphValuesWorkbook);
reset_output_file(prismWorkbook);

cbPositivePairedTable = pairedTable(pairedTable.CbStatus == "CB+", :);

write_table_sheet(graphValuesWorkbook, 'README', make_export_readme_table( ...
    'Verification workbook', data, opts, ...
    'Each numbered sheet contains the mouse-level values, counts, and denominators used for the matching figure panel.', ...
    layerAllocationSource));
write_table_sheet(graphValuesWorkbook, '01_layer_allocation', layerAllocationTable);
write_table_sheet(graphValuesWorkbook, '02_status_overall', statusOverallTable);
write_table_sheet(graphValuesWorkbook, '03_status_composite', statusCompositeTable);
write_table_sheet(graphValuesWorkbook, '04_cb_positive_layer', cbPositivePairedTable);
write_table_sheet(graphValuesWorkbook, '05_status_split', pairedTable);
write_table_sheet(graphValuesWorkbook, '06_stacked_composition', summaryTable);
write_table_sheet(graphValuesWorkbook, '07_heatmap_within', layerCompositeWithinTable);
write_table_sheet(graphValuesWorkbook, '08_heatmap_global', layerCompositeGlobalTable);

write_table_sheet(prismWorkbook, 'README', make_export_readme_table( ...
    'Prism-ready workbook', data, opts, ...
    'Sheets are reshaped into wide tables for import or copy/paste into GraphPad Prism.', ...
    layerAllocationSource));
write_table_sheet(prismWorkbook, 'P1_layer_allocation', make_layer_allocation_prism_table(layerAllocationTable));
write_table_sheet(prismWorkbook, 'P2_status_overall', make_status_overall_prism_table(statusOverallTable));
write_table_sheet(prismWorkbook, 'P3_composite_raw', make_status_composite_prism_table(statusCompositeTable));
write_table_sheet(prismWorkbook, 'P3_composite_stack', make_status_composite_stack_prism_table(statusCompositeTable));
write_table_sheet(prismWorkbook, 'P4_cb_positive_layer', make_cb_positive_layer_prism_table(pairedTable));
write_table_sheet(prismWorkbook, 'P5_status_split', make_status_split_prism_table(pairedTable));
write_table_sheet(prismWorkbook, 'P6_stacked_composition', make_stacked_composition_prism_table(summaryTable));
write_table_sheet(prismWorkbook, 'P7_heatmap_within', make_composite_heatmap_prism_table(layerCompositeWithinTable));
write_table_sheet(prismWorkbook, 'P8_heatmap_global', make_composite_heatmap_prism_table(layerCompositeGlobalTable));
end

function readmeTable = make_export_readme_table(workbookKind, data, opts, description, layerAllocationSource)
Item = {
    'Workbook type'
    'Description'
    'Source workbook'
    'Source sheet'
    'Panel 1 source sheet'
    'Side filter'
    'Generated'
    'Denominator'
    'Limitation'};
Value = {
    workbookKind
    description
    char(data.workbookFile)
    char(data.sheetName)
    char(layerAllocationSource)
    char(opts.SideFilter)
    datestr(now, 31)
    'mKate/RFP-tagged cells unless the sheet or panel title states otherwise'
    'These data do not contain total CB+ or CB- population denominators independent of mKate/RFP tagging'};
readmeTable = table(Item, Value, 'VariableNames', {'Item', 'Value'});
end

function write_table_sheet(workbookFile, sheetName, sheetTable)
writetable(sheetTable, workbookFile, 'Sheet', sheetName);
end

function reset_output_file(fileName)
if exist(fileName, 'file') == 2
    delete(fileName);
end
end

function prismTable = make_layer_allocation_prism_table(layerAllocationTable)
Mouse = layerAllocationTable.Mouse;
Superficial = layerAllocationTable.SuperficialPercentOfInLayer;
Deep = layerAllocationTable.DeepPercentOfInLayer;
prismTable = table(Mouse, Superficial, Deep);
end

function prismTable = make_status_overall_prism_table(statusOverallTable)
mice = unique(statusOverallTable.Mouse);
Mouse = mice(:);
CBPositive = collect_status_overall_values(statusOverallTable, mice, "CB+");
CBNegative = collect_status_overall_values(statusOverallTable, mice, "CB-");
CBIntermediate = collect_status_overall_values(statusOverallTable, mice, "CB+/-");
prismTable = table(Mouse, CBPositive, CBNegative, CBIntermediate);
end

function values = collect_status_overall_values(statusOverallTable, mice, statusName)
values = nan(numel(mice), 1);
for mouseIdx = 1:numel(mice)
    row = statusOverallTable.Mouse == mice(mouseIdx) & statusOverallTable.CbStatus == statusName;
    values(mouseIdx) = first_value(statusOverallTable.PercentOfInLayerMKate(row));
end
end

function prismTable = make_status_composite_prism_table(statusCompositeTable)
Mouse = statusCompositeTable.Mouse;
CBNegative = statusCompositeTable.CbNegativePercent;
CBPositiveIntermediate = statusCompositeTable.CbPositiveIntermediatePercent;
CBPositive = statusCompositeTable.CbPositivePercent;
CBIntermediate = statusCompositeTable.CbIntermediatePercent;
prismTable = table(Mouse, CBNegative, CBPositiveIntermediate, CBPositive, CBIntermediate);
end

function prismTable = make_status_composite_stack_prism_table(statusCompositeTable)
Component = {'CB-'; 'CB+'; 'CB+/- intermediate'};
CBNegativeColumn = [mean_omit_nan(statusCompositeTable.CbNegativePercent); 0; 0];
CBPositiveIntermediateColumn = [0; ...
    mean_omit_nan(statusCompositeTable.CbPositivePercent); ...
    mean_omit_nan(statusCompositeTable.CbIntermediatePercent)];
prismTable = table(Component, CBNegativeColumn, CBPositiveIntermediateColumn);
end

function prismTable = make_cb_positive_layer_prism_table(pairedTable)
cbPositiveTable = pairedTable(pairedTable.CbStatus == "CB+", :);
Mouse = cbPositiveTable.Mouse;
Superficial = cbPositiveTable.SuperficialPercent;
Deep = cbPositiveTable.DeepPercent;
prismTable = table(Mouse, Superficial, Deep);
end

function prismTable = make_status_split_prism_table(pairedTable)
mice = unique(pairedTable.Mouse);
Mouse = mice(:);
CBPositiveSuperficial = collect_paired_status_values(pairedTable, mice, "CB+", 'SuperficialPercent');
CBPositiveDeep = collect_paired_status_values(pairedTable, mice, "CB+", 'DeepPercent');
CBNegativeSuperficial = collect_paired_status_values(pairedTable, mice, "CB-", 'SuperficialPercent');
CBNegativeDeep = collect_paired_status_values(pairedTable, mice, "CB-", 'DeepPercent');
CBIntermediateSuperficial = collect_paired_status_values(pairedTable, mice, "CB+/-", 'SuperficialPercent');
CBIntermediateDeep = collect_paired_status_values(pairedTable, mice, "CB+/-", 'DeepPercent');
prismTable = table(Mouse, CBPositiveSuperficial, CBPositiveDeep, CBNegativeSuperficial, ...
    CBNegativeDeep, CBIntermediateSuperficial, CBIntermediateDeep);
end

function values = collect_paired_status_values(pairedTable, mice, statusName, columnName)
values = nan(numel(mice), 1);
for mouseIdx = 1:numel(mice)
    row = pairedTable.Mouse == mice(mouseIdx) & pairedTable.CbStatus == statusName;
    values(mouseIdx) = first_value(pairedTable.(columnName)(row));
end
end

function prismTable = make_stacked_composition_prism_table(summaryTable)
mice = unique(summaryTable.Mouse);
layers = ["Superficial"; "Deep"];
nRows = numel(mice) * numel(layers);
Mouse = nan(nRows, 1);
Layer = strings(nRows, 1);
CBPositive = nan(nRows, 1);
CBNegative = nan(nRows, 1);
CBIntermediate = nan(nRows, 1);
LayerTotal = nan(nRows, 1);

k = 0;
for mouseIdx = 1:numel(mice)
    for layerIdx = 1:numel(layers)
        k = k + 1;
        mouse = mice(mouseIdx);
        layer = layers(layerIdx);
        Mouse(k) = mouse;
        Layer(k) = layer;
        CBPositive(k) = get_layer_status_percent(summaryTable, mouse, layer, "CB+");
        CBNegative(k) = get_layer_status_percent(summaryTable, mouse, layer, "CB-");
        CBIntermediate(k) = get_layer_status_percent(summaryTable, mouse, layer, "CB+/-");
        LayerTotal(k) = get_layer_total_count(summaryTable, mouse, layer);
    end
end

prismTable = table(Mouse, Layer, CBPositive, CBNegative, CBIntermediate, LayerTotal);
end

function value = get_layer_status_percent(summaryTable, mouse, layer, cbStatus)
row = summaryTable.Mouse == mouse & summaryTable.Layer == layer & summaryTable.CbStatus == cbStatus;
value = first_value(summaryTable.PercentOfLayer(row));
end

function value = get_layer_total_count(summaryTable, mouse, layer)
row = summaryTable.Mouse == mouse & summaryTable.Layer == layer;
value = first_value(summaryTable.TotalLayerCount(row));
end

function prismTable = make_composite_heatmap_prism_table(compositeTable)
[heatValues, mice] = make_composite_heatmap_matrix(compositeTable);
Mouse = mice(:);
Superficial_CBPositiveIntermediate = heatValues(:, 1);
Superficial_CBNegative = heatValues(:, 2);
Deep_CBPositiveIntermediate = heatValues(:, 3);
Deep_CBNegative = heatValues(:, 4);
prismTable = table(Mouse, Superficial_CBPositiveIntermediate, Superficial_CBNegative, ...
    Deep_CBPositiveIntermediate, Deep_CBNegative);
end

function write_panel_guide(outputFolder, data, opts, layerAllocationSource, layerAllocationMouseCount)
guideFile = fullfile(outputFolder, 'cb_quantification_panel_guide.txt');
fileId = fopen(guideFile, 'w');
if fileId < 0
    warning('CBPlot:PanelGuideWriteFailed', 'Could not write panel guide: %s', guideFile);
    return
end

rawTable = data.raw;
mouseCount = numel(unique(rawTable.Mouse(~isnan(rawTable.Mouse))));
stackCount = height(rawTable);

fprintf(fileId, 'CB image quantification panel guide\n');
fprintf(fileId, '=================================\n\n');
fprintf(fileId, 'Biological framing\n');
fprintf(fileId, 'Recent evidence suggests SWRs differentiate deep and superficial CA1. These panels ask whether the mKate/RFP-tagged engram is allocated differently across deep versus superficial CA1, and whether tagged cells are preferentially CB+ versus CB-.\n\n');
fprintf(fileId, 'Input and denominator\n');
fprintf(fileId, 'Input workbook: %s\n', char(data.workbookFile));
fprintf(fileId, 'CB-status sheet: %s\n', char(data.sheetName));
fprintf(fileId, 'Panel 1 source sheet: %s\n', char(layerAllocationSource));
fprintf(fileId, 'Side filter: %s\n', char(opts.SideFilter));
fprintf(fileId, 'Imported stacks: %d\n', stackCount);
fprintf(fileId, 'CB-status mice: %d\n', mouseCount);
fprintf(fileId, 'Panel 1 mice: %d\n', layerAllocationMouseCount);
fprintf(fileId, 'All percentages use mKate/RFP-tagged cells as the denominator unless noted otherwise.\n\n');
fprintf(fileId, 'Excel exports\n');
fprintf(fileId, 'cb_graph_values_for_verification.xlsx: exact mouse-level values, counts, and denominators used for each graph.\n');
fprintf(fileId, 'cb_prism_ready_tables.xlsx: wide Prism-ready tables for importing or copying into GraphPad Prism.\n\n');
fprintf(fileId, 'Panels\n');
fprintf(fileId, '1. cb_mkate_layer_allocation_paired: %% of in-layer mKate/RFP cells in superficial versus deep CA1. This uses %s and addresses whether the tagged engram is preferentially allocated to one CA1 layer.\n', char(layerAllocationSource));
fprintf(fileId, '2. cb_status_overall_assignment: %% of in-layer mKate/RFP cells assigned to CB+, CB-, or CB+/- categories after pooling superficial and deep CA1. This addresses CB+ versus CB- assignment overall.\n');
fprintf(fileId, '3. cb_status_composite_assignment_column: column plot contrasting CB- with a stacked CB+ plus CB+/- intermediate composite. Individual mouse points show CB- versus composite totals.\n');
fprintf(fileId, '4. cb_positive_deep_vs_superficial_paired: per-layer %% of mKate/RFP cells that are CB+. This addresses whether CB+ assignment differs between superficial and deep CA1.\n');
fprintf(fileId, '5. cb_status_deep_vs_superficial_paired: per-layer %% of mKate/RFP cells in each CB status category.\n');
fprintf(fileId, '6. cb_status_stacked_composition_by_mouse_layer: stacked composition of CB status within each mouse/layer. Each bar sums to 100%% of mKate/RFP cells in that mouse/layer.\n');
fprintf(fileId, '7. cb_status_layer_percent_heatmap: rows are mice; columns are superficial or deep crossed with CB+/intermediate or CB-. Tile color and text show %% of mKate/RFP cells within the same mouse and layer. For each mouse, S CB+/int + S CB- sums to 100%%, and D CB+/int + D CB- sums to 100%%.\n');
fprintf(fileId, '8. cb_status_layer_global_percent_heatmap: rows are mice; columns are superficial or deep crossed with CB+/intermediate or CB-. Here all four columns sum to 100%% for each mouse.\n\n');
fprintf(fileId, 'Current limitation\n');
fprintf(fileId, 'These data cannot answer: out of all CB+ cells or all CB- cells, how many are mKate/RFP tagged? That would require total CB+ and CB- population denominators independent of mKate/RFP tagging.\n');

fclose(fileId);
end

function write_overview_html(outputFolder, data, opts, layerAllocationSource, layerAllocationMouseCount)
overviewFile = fullfile(outputFolder, 'overview.html');
fileId = fopen(overviewFile, 'w');
if fileId < 0
    warning('CBPlot:OverviewWriteFailed', 'Could not write overview: %s', overviewFile);
    return
end

rawTable = data.raw;
mouseCount = numel(unique(rawTable.Mouse(~isnan(rawTable.Mouse))));
stackCount = height(rawTable);

fprintf(fileId, '<!doctype html>\n');
fprintf(fileId, '<html lang="en">\n<head>\n<meta charset="utf-8">\n');
fprintf(fileId, '<meta name="viewport" content="width=device-width, initial-scale=1">\n');
fprintf(fileId, '<title>CB Image Quantification Overview</title>\n');
fprintf(fileId, '<style>\n');
fprintf(fileId, 'body{margin:0;font-family:Arial,Helvetica,sans-serif;background:#f6f7f8;color:#202124;line-height:1.45;}\n');
fprintf(fileId, 'header{background:#ffffff;border-bottom:1px solid #d9dde3;padding:28px 36px;}\n');
fprintf(fileId, 'main{max-width:1180px;margin:0 auto;padding:28px 24px 48px;}\n');
fprintf(fileId, 'h1{margin:0 0 10px;font-size:28px;font-weight:700;}\n');
fprintf(fileId, 'h2{margin:0 0 10px;font-size:20px;}\n');
fprintf(fileId, 'p{margin:0 0 10px;}\n');
fprintf(fileId, '.meta{color:#56606b;font-size:14px;}\n');
fprintf(fileId, '.note{background:#fff7e6;border:1px solid #f0d59a;border-radius:6px;padding:14px 16px;margin:18px 0;}\n');
fprintf(fileId, '.panel{background:#ffffff;border:1px solid #d9dde3;border-radius:6px;margin:0 0 26px;overflow:hidden;}\n');
fprintf(fileId, '.panel-body{padding:18px 20px 20px;}\n');
fprintf(fileId, '.figure-wrap{background:#ffffff;border-bottom:1px solid #e5e8ec;text-align:center;padding:16px;}\n');
fprintf(fileId, 'img{max-width:100%%;height:auto;border:1px solid #edf0f2;}\n');
fprintf(fileId, '.links{margin-top:10px;font-size:14px;}\n');
fprintf(fileId, '.links a{color:#1b5e9e;text-decoration:none;margin-right:14px;}\n');
fprintf(fileId, '.label{font-weight:700;color:#374151;}\n');
fprintf(fileId, 'footer{max-width:1180px;margin:0 auto;padding:0 24px 36px;color:#667085;font-size:13px;}\n');
fprintf(fileId, '</style>\n</head>\n<body>\n');

fprintf(fileId, '<header>\n');
fprintf(fileId, '<h1>CB Image Quantification Overview</h1>\n');
fprintf(fileId, '<p>These panels ask whether the mKate/RFP-tagged engram is allocated differently across superficial versus deep CA1, and whether tagged cells are preferentially assigned to CB+ versus CB- categories.</p>\n');
fprintf(fileId, '<p class="meta">Workbook: %s<br>CB-status sheet: %s<br>Panel 1 source sheet: %s<br>Side filter: %s<br>Imported stacks: %d, CB-status mice: %d, panel 1 mice: %d<br>Generated: %s</p>\n', ...
    html_escape(data.workbookFile), html_escape(data.sheetName), html_escape(layerAllocationSource), html_escape(opts.SideFilter), ...
    stackCount, mouseCount, layerAllocationMouseCount, datestr(now, 31));
fprintf(fileId, '</header>\n<main>\n');

fprintf(fileId, '<div class="note"><p><span class="label">Denominator:</span> All percentages use mKate/RFP-tagged cells as the denominator unless a panel says otherwise.</p>');
fprintf(fileId, '<p><span class="label">Limitation:</span> These data cannot answer what fraction of all CB+ or all CB- cells are mKate/RFP-tagged, because the sheet does not contain total CB+ and CB- population denominators independent of mKate/RFP tagging.</p>');
fprintf(fileId, '<p class="links"><span class="label">Data exports:</span> <a href="cb_graph_values_for_verification.xlsx">Verification workbook</a><a href="cb_prism_ready_tables.xlsx">Prism-ready workbook</a></p></div>\n');

write_overview_panel(fileId, 'cb_mkate_layer_allocation_paired', ...
    'Engram allocation across CA1 layers', ...
    'Question: Are mKate/RFP-tagged cells distributed differently between superficial and deep CA1?', ...
    sprintf('Interpretation: Each connected pair is one mouse, using %s. The two y-values are the percent of all in-layer mKate/RFP cells found in superficial versus deep CA1. A consistent shift in one direction suggests layer-biased allocation of the tagged engram.', char(layerAllocationSource)));

write_overview_panel(fileId, 'cb_status_overall_assignment', ...
    'Overall assignment to CB status', ...
    'Question: Are mKate/RFP-tagged cells preferentially assigned to CB+ rather than CB- categories?', ...
    'Interpretation: Values are pooled across superficial and deep CA1 for each mouse. Higher CB+ values indicate more mKate/RFP-tagged cells are CB+; higher CB- values indicate the tagged population is mostly CB-. CB+/- is kept separate as an ambiguous/mixed category.');

write_overview_panel(fileId, 'cb_status_composite_assignment_column', ...
    'Composite CB assignment column plot', ...
    'Question: If CB+ and CB+/- intermediate cells are considered together, how does that composite compare with CB- assignment?', ...
    'Interpretation: The CB- column is a single category. The composite column is stacked to show the CB+ contribution plus the CB+/- intermediate contribution. Individual mouse points compare CB- against the combined CB+/intermediate total.');

write_overview_panel(fileId, 'cb_positive_deep_vs_superficial_paired', ...
    'CB+ fraction by CA1 layer', ...
    'Question: Within each CA1 layer, what percent of mKate/RFP-tagged cells are CB+?', ...
    'Interpretation: Each mouse contributes a superficial and deep value. A lower deep value than superficial value means CB+ assignment among tagged cells is lower in deep CA1 for that mouse.');

write_overview_panel(fileId, 'cb_status_deep_vs_superficial_paired', ...
    'Per-layer assignment to CB status', ...
    'Question: How do CB+, CB-, and CB+/- fractions differ between superficial and deep CA1?', ...
    'Interpretation: This expands the CB+ panel to all status categories. Read each small plot separately: each mouse has paired superficial and deep percentages for that CB status.');

write_overview_panel(fileId, 'cb_status_stacked_composition_by_mouse_layer', ...
    'CB status composition by mouse and layer', ...
    'Question: What is the full CB-status composition of tagged cells in each mouse/layer?', ...
    'Interpretation: Each bar is one mouse-layer combination and sums to 100% of mKate/RFP-tagged cells in that layer. This is useful for seeing whether CB- dominates deep CA1 or whether a mouse has an unusual CB+/- fraction.');

write_overview_panel(fileId, 'cb_status_layer_percent_heatmap', ...
    'Within-layer composite heatmap', ...
    'Question: Which mice drive the layer-specific CB-status patterns?', ...
    'Interpretation: Rows are mice. Columns are superficial or deep crossed with CB+/intermediate or CB-. Tile color and text show percent of mKate/RFP-tagged cells within the same mouse and layer. For each mouse, the two superficial columns sum to 100%, and the two deep columns sum to 100%.');

write_overview_panel(fileId, 'cb_status_layer_global_percent_heatmap', ...
    'Global layer/status composite heatmap', ...
    'Question: How are tagged cells distributed across both layer and CB status when all superficial and deep mKate/RFP cells are the denominator?', ...
    'Interpretation: Rows are mice. The four columns are superficial CB+/intermediate, superficial CB-, deep CB+/intermediate, and deep CB-. For each mouse, all four columns sum to 100%, so this shows layer allocation and CB assignment in one denominator.');

fprintf(fileId, '</main>\n');
fprintf(fileId, '<footer>Editable vector PDFs, source CSV tables, and Excel export workbooks are saved in the same folder as this overview.</footer>\n');
fprintf(fileId, '</body>\n</html>\n');
fclose(fileId);
end

function write_overview_panel(fileId, baseName, titleText, questionText, interpretationText)
pngName = [baseName '.png'];
pdfName = [baseName '.pdf'];
fprintf(fileId, '<section class="panel">\n');
fprintf(fileId, '<div class="figure-wrap"><img src="%s" alt="%s"></div>\n', ...
    html_escape(pngName), html_escape(titleText));
fprintf(fileId, '<div class="panel-body">\n');
fprintf(fileId, '<h2>%s</h2>\n', html_escape(titleText));
fprintf(fileId, '<p><span class="label">%s</span></p>\n', html_escape(questionText));
fprintf(fileId, '<p>%s</p>\n', html_escape(interpretationText));
fprintf(fileId, '<p class="links"><a href="%s">PNG</a><a href="%s">Editable PDF</a></p>\n', ...
    html_escape(pngName), html_escape(pdfName));
fprintf(fileId, '</div>\n</section>\n');
end

function out = html_escape(value)
out = char(string(value));
out = strrep(out, '&', '&amp;');
out = strrep(out, '<', '&lt;');
out = strrep(out, '>', '&gt;');
out = strrep(out, '"', '&quot;');
out = strrep(out, '''', '&#39;');
end

function colors = cb_status_colors()
colors = [
    0.74 0.19 0.22
    0.18 0.37 0.67
    0.46 0.46 0.46];
end

function label = title_with_filter(baseTitle, opts)
if lower(opts.SideFilter) == "all"
    label = baseTitle;
else
    label = sprintf('%s | %s', baseTitle, char(opts.SideFilter));
end
end

function value = first_value(values)
if isempty(values)
    value = nan;
else
    value = values(1);
end
end

function value = safe_percent(numerator, denominator)
if denominator == 0
    value = nan;
else
    value = numerator * 100 / denominator;
end
end

function shortName = layer_short_name(layerName)
if layerName == "Superficial"
    shortName = 'S';
elseif layerName == "Deep"
    shortName = 'D';
else
    shortName = char(layerName);
end
end

function outputFolder = default_output_folder()
codeFolder = fileparts(mfilename('fullpath'));
repoRoot = fullfile(codeFolder, '..', '..');
candidateFolders = {
    fullfile(repoRoot, 'Results')
    fullfile(repoRoot, '..', 'Results')};
outputFolder = fullfile(first_existing_folder(candidateFolders), 'CB_quantification_figures');
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
validNames = {'outputfolder', 'sidefilter', 'savefigures', 'closefigures', ...
    'cachefile', 'reimportexcel'};
tf = any(strcmp(name, validNames));
end
