function data = import_interneuron_quantifications(workbookFile, varargin)
%IMPORT_INTERNEURON_QUANTIFICATIONS Load interneuron marker quantifications.
%   DATA = IMPORT_INTERNEURON_QUANTIFICATIONS(WORKBOOKFILE) imports the
%   interneuron sheets from Quantifications.xlsx and returns raw count
%   records, active/total measurements, animal-level summaries, and
%   denominator availability tables.
%
%   Imported datasets are kept separate:
%     intrinsic_cFos    cFos immunostaining in endogenous cFos experiments.
%     cfos_TRE_mKate    mKate/RFP-tagged cells in cfos-TRE + AAV-TTA mice.
%
%   Name-value inputs:
%     IncludeZeroCounts  Keep zero-count active rows when true.
%     CacheFile          MAT cache path for imported workbook data.
%     ReimportExcel      Ignore the cache and re-read the workbook.
%
%   Count summaries include raw counts and normalized count columns. The
%   normalized count columns use 60x-equivalent stacks: 20x mouse-level
%   aSNCG+cFos+CCK8 rows 3:6 and 17:20 use 6 stacks and a 3x FOV factor.
%   RFP/mKate denominators use explicit RFP-negative blocks when present.
%   cfos_TRE_mKate slice groups with no mKate/RFP-positive cells or no
%   positive interneuron-marker cells are excluded before summaries.

if nargin < 1
    workbookFile = [];
elseif is_name_value_start(workbookFile)
    varargin = [{workbookFile}, varargin];
    workbookFile = [];
end

if isempty(workbookFile)
    workbookFile = default_workbook_file();
end

parser = inputParser;
addParameter(parser, 'IncludeZeroCounts', true, @(x) islogical(x) && isscalar(x));
addParameter(parser, 'CacheFile', '', @(x) ischar(x) || isstring(x));
addParameter(parser, 'ReimportExcel', false, @(x) islogical(x) && isscalar(x));
parse(parser, varargin{:});

includeZeroCounts = parser.Results.IncludeZeroCounts;
cacheFile = char(parser.Results.CacheFile);
reimportExcel = parser.Results.ReimportExcel;

if isempty(cacheFile)
    cacheFile = default_cache_file(workbookFile);
end

if ~reimportExcel && isfile(cacheFile)
    loaded = load(cacheFile, 'data');
    data = loaded.data;
    if cache_has_required_fields(data) || ~isfile(workbookFile)
        data = ensure_exclusion_fields(data);
        data.loadedFrom = "mat";
        data = apply_include_zero_counts(data, includeZeroCounts);
        return
    end
end

data = import_from_excel(workbookFile);
data.cacheFile = string(cacheFile);
data.loadedFrom = "excel";
data.importedAt = string(datetime('now'));

cacheFolder = fileparts(cacheFile);
if ~isempty(cacheFolder) && ~exist(cacheFolder, 'dir')
    mkdir(cacheFolder);
end
save(cacheFile, 'data');

data = apply_include_zero_counts(data, includeZeroCounts);
end

function data = import_from_excel(workbookFile)
if ~isfile(workbookFile)
    error('InterneuronImport:MissingWorkbook', 'Workbook not found: %s', char(workbookFile));
end

availableSheets = string(sheetnames(workbookFile));
configs = block_configs();
records = empty_record_table();

for cfgIdx = 1:numel(configs)
    cfg = configs(cfgIdx);
    if ~any(strcmpi(availableSheets, cfg.Sheet))
        warning('InterneuronImport:MissingSheet', ...
            'Sheet "%s" not found in %s.', cfg.Sheet, char(workbookFile));
        continue
    end

    rawCells = readcell(workbookFile, 'Sheet', cfg.Sheet, 'Range', cfg.Range);
    if cfg.Mode == "regionBlock"
        records = [records; import_region_block(rawCells, cfg)]; %#ok<AGROW>
    elseif cfg.Mode == "allRegionTotals"
        records = [records; import_all_region_totals(rawCells, cfg)]; %#ok<AGROW>
    else
        error('InterneuronImport:UnknownMode', 'Unknown import mode "%s".', cfg.Mode);
    end
end

records = normalize_record_types(records);
[records, excludedSlices, excludedRows] = exclude_empty_mkate_slices(records);
measurements = make_measurements(records);

data = struct();
data.importerVersion = 4;
data.workbookFile = string(workbookFile);
data.records = records;
data.excludedSlices = excludedSlices;
data.excludedRows = excludedRows;
data.measurements = measurements;
data.byAnimalRegion = summarize_measurements(measurements, ...
    {'Dataset', 'MouseLine', 'Reporter', 'Sheet', 'MarkerPanel', 'Mouse', 'Marker', 'Region'});
data.byAnimalMarker = summarize_measurements(measurements, ...
    {'Dataset', 'MouseLine', 'Reporter', 'Sheet', 'MarkerPanel', 'Mouse', 'Marker'});
data.byDatasetMarkerRegion = summarize_measurements(measurements, ...
    {'Dataset', 'MouseLine', 'Reporter', 'Sheet', 'MarkerPanel', 'Marker', 'Region'});
data.totalAvailability = summarize_total_availability(measurements);
data.regionOrder = ["S Oriens"; "Deep S Pyr"; "Superficial S Pyr"; "S Radiatum"];
data.notes = make_notes();
end

function configs = block_configs()
regions = ["S Oriens", "S Radiatum", "Deep S Pyr", "Superficial S Pyr"];

configs = repmat(empty_config(), 0, 1);

configs(end + 1) = make_region_config( ...
    "PV+SOM+RFP", "PV/SOM", "cfos_TRE_mKate", "cfos-TRE + AAV-TTA", "mKate/RFP", ...
    "Active", "region", "A:CR", 3, 301, "B", "C", "A", "D", ...
    ["E", "I", "M", "Q"], 4, regions, ...
    ["PV", "SST/SOM", "PV+SST/SOM", "PV-/SST-"]);

configs(end + 1) = make_region_config( ...
    "PV+SOM+RFP", "PV/SOM", "cfos_TRE_mKate", "cfos-TRE + AAV-TTA", "mKate/RFP", ...
    "Total", "region", "A:CR", 3, 301, "BK", "BL", "BJ", "BM", ...
    ["BN", "BQ", "BT", "BW"], 3, regions, ...
    ["PV", "SST/SOM", "PV+SST/SOM"]);

configs(end + 1) = make_region_config( ...
    "PV+SOM+RFP", "PV/SOM", "cfos_TRE_mKate", "cfos-TRE + AAV-TTA", "mKate/RFP", ...
    "Inactive", "region", "A:DO", 3, 301, "CP", "CQ", "CO", "CR", ...
    ["CS", "CV", "CY", "DB"], 3, regions, ...
    ["PV", "SST/SOM", "PV+SST/SOM"]);

configs(end + 1) = make_region_config( ...
    "PV+CFOS+SOM", "PV/SOM", "intrinsic_cFos", "cFos immunostaining", "cFos", ...
    "Active", "all_regions", "A:AS", 3, 48, "A", "B", "", "C", ...
    ["D", "H", "L", "P"], 4, regions, ...
    ["PV", "SST/SOM", "PV+SST/SOM", "PV-/SST-"]);

configs(end + 1) = make_region_config( ...
    "PV+CFOS+SOM", "PV/SOM", "intrinsic_cFos", "cFos immunostaining", "cFos", ...
    "Inactive", "region", "A:BJ", 3, 42, "AW", "AX", "", "", ...
    ["AY", "BB", "BE", "BH"], 3, regions, ...
    ["PV", "SST/SOM", "PV+SST/SOM"]);

configs(end + 1) = make_total_config( ...
    "PV+CFOS+SOM", "PV/SOM", "intrinsic_cFos", "cFos immunostaining", "cFos", ...
    "A:AS", 11, 14, "AP", ["AQ", "AR", "AS"], ...
    ["PV", "SST/SOM", "PV+SST/SOM"]);

configs(end + 1) = make_region_config( ...
    "aSNCG+cFos+CCK8", "SNCG/CCK8", "intrinsic_cFos", "cFos immunostaining", "cFos", ...
    "Active", "all_regions", "A:CQ", 3, 6, "A", "", "", "", ...
    ["E", "I", "M", "Q"], 4, regions, ...
    ["SNCG", "CCK8", "SNCG+CCK8", "SNCG-/CCK8-"], 6, 3);

configs(end + 1) = make_region_config( ...
    "aSNCG+cFos+CCK8", "SNCG/CCK8", "intrinsic_cFos", "cFos immunostaining", "cFos", ...
    "Inactive", "region", "A:CQ", 3, 6, "AO", "", "", "", ...
    ["AP", "AS", "AV", "AY"], 3, regions, ...
    ["SNCG", "CCK8", "SNCG+CCK8"], 6, 3);

configs(end + 1) = make_total_config( ...
    "aSNCG+cFos+CCK8", "SNCG/CCK8", "intrinsic_cFos", "cFos immunostaining", "cFos", ...
    "A:CQ", 3, 6, "AC", ["AD", "AE", "AF"], ...
    ["SNCG", "CCK8", "SNCG+CCK8"], 6, 3);

configs(end + 1) = make_region_config( ...
    "aSNCG+cFos+CCK8", "proCCK/VIP", "intrinsic_cFos", "cFos immunostaining", "cFos", ...
    "Active", "all_regions", "A:CQ", 17, 20, "A", "", "", "", ...
    ["E", "I", "M", "Q"], 4, regions, ...
    ["proCCK", "VIP", "proCCK+VIP", "proCCK-/VIP-"], 6, 3);

configs(end + 1) = make_region_config( ...
    "aSNCG+cFos+CCK8", "proCCK/VIP", "intrinsic_cFos", "cFos immunostaining", "cFos", ...
    "Inactive", "region", "A:CQ", 17, 20, "AO", "", "", "", ...
    ["AP", "AS", "AV", "AY"], 3, regions, ...
    ["proCCK", "VIP", "proCCK+VIP"], 6, 3);

configs(end + 1) = make_total_config( ...
    "aSNCG+cFos+CCK8", "proCCK/VIP", "intrinsic_cFos", "cFos immunostaining", "cFos", ...
    "A:CQ", 17, 20, "AC", ["AD", "AE", "AF"], ...
    ["proCCK", "VIP", "proCCK+VIP"], 6, 3);

configs(end + 1) = make_region_config( ...
    "aSNCG+cFos+CCK8", "proCCK/VIP", "cfos_TRE_mKate", "cfos-TRE + AAV-TTA", "mKate/RFP", ...
    "Active", "region", "A:DX", 38, 279, "B", "C", "A", "", ...
    ["E", "I", "M", "Q"], 4, regions, ...
    ["proCCK", "VIP", "proCCK+VIP", "proCCK-/VIP-"]);

configs(end + 1) = make_region_config( ...
    "aSNCG+cFos+CCK8", "proCCK/VIP", "cfos_TRE_mKate", "cfos-TRE + AAV-TTA", "mKate/RFP", ...
    "Total", "region", "A:DX", 38, 279, "BJ", "BK", "BI", "", ...
    ["BM", "BQ", "BU", "BY"], 3, regions, ...
    ["proCCK", "VIP", "proCCK+VIP"]);

configs(end + 1) = make_region_config( ...
    "aSNCG+cFos+CCK8", "proCCK/VIP", "cfos_TRE_mKate", "cfos-TRE + AAV-TTA", "mKate/RFP", ...
    "Inactive", "region", "A:DS", 42, 283, "CZ", "DA", "CY", "", ...
    ["DC", "DF", "DI", "DL"], 3, regions, ...
    ["proCCK", "VIP", "proCCK+VIP"]);
end

function cfg = empty_config()
cfg = struct();
cfg.Mode = "";
cfg.Sheet = "";
cfg.MarkerPanel = "";
cfg.Dataset = "";
cfg.MouseLine = "";
cfg.Reporter = "";
cfg.CountRole = "";
cfg.DenominatorScope = "";
cfg.Range = "";
cfg.FirstRow = nan;
cfg.LastRow = nan;
cfg.MouseCol = "";
cfg.StackCol = "";
cfg.SideCol = "";
cfg.AddressCol = "";
cfg.LayerStartCols = strings(0, 1);
cfg.MarkersPerLayer = nan;
cfg.Regions = strings(0, 1);
cfg.Markers = strings(0, 1);
cfg.TotalCols = strings(0, 1);
cfg.SlicesRepresented = 1;
cfg.FovAreaFactor = 1;
end

function cfg = make_region_config(sheet, markerPanel, dataset, mouseLine, reporter, ...
    countRole, denominatorScope, rangeName, firstRow, lastRow, mouseCol, stackCol, sideCol, ...
    addressCol, layerStartCols, markersPerLayer, regions, markers, slicesRepresented, fovAreaFactor)
if nargin < 19 || isempty(slicesRepresented)
    slicesRepresented = 1;
end
if nargin < 20 || isempty(fovAreaFactor)
    fovAreaFactor = 1;
end
cfg = empty_config();
cfg.Mode = "regionBlock";
cfg.Sheet = sheet;
cfg.MarkerPanel = markerPanel;
cfg.Dataset = dataset;
cfg.MouseLine = mouseLine;
cfg.Reporter = reporter;
cfg.CountRole = countRole;
cfg.DenominatorScope = denominatorScope;
cfg.Range = rangeName;
cfg.FirstRow = firstRow;
cfg.LastRow = lastRow;
cfg.MouseCol = mouseCol;
cfg.StackCol = stackCol;
cfg.SideCol = sideCol;
cfg.AddressCol = addressCol;
cfg.LayerStartCols = layerStartCols(:);
cfg.MarkersPerLayer = markersPerLayer;
cfg.Regions = regions(:);
cfg.Markers = markers(:);
cfg.SlicesRepresented = slicesRepresented;
cfg.FovAreaFactor = fovAreaFactor;
end

function cfg = make_total_config(sheet, markerPanel, dataset, mouseLine, reporter, ...
    rangeName, firstRow, lastRow, mouseCol, totalCols, markers, slicesRepresented, fovAreaFactor)
if nargin < 12 || isempty(slicesRepresented)
    slicesRepresented = 1;
end
if nargin < 13 || isempty(fovAreaFactor)
    fovAreaFactor = 1;
end
cfg = empty_config();
cfg.Mode = "allRegionTotals";
cfg.Sheet = sheet;
cfg.MarkerPanel = markerPanel;
cfg.Dataset = dataset;
cfg.MouseLine = mouseLine;
cfg.Reporter = reporter;
cfg.CountRole = "Total";
cfg.DenominatorScope = "all_regions";
cfg.Range = rangeName;
cfg.FirstRow = firstRow;
cfg.LastRow = lastRow;
cfg.MouseCol = mouseCol;
cfg.TotalCols = totalCols(:);
cfg.Markers = markers(:);
cfg.SlicesRepresented = slicesRepresented;
cfg.FovAreaFactor = fovAreaFactor;
end

function records = import_region_block(rawCells, cfg)
records = empty_record_table();
mouseCol = excel_col_to_num(cfg.MouseCol);
stackCol = excel_col_to_num(cfg.StackCol);
sideCol = excel_col_to_num(cfg.SideCol);
addressCol = excel_col_to_num(cfg.AddressCol);
layerStartCols = arrayfun(@excel_col_to_num, cfg.LayerStartCols);
lastSide = "";

for rowIdx = cfg.FirstRow:min(cfg.LastRow, size(rawCells, 1))
    mouse = cell_to_double(value_at(rawCells, rowIdx, mouseCol));
    if isnan(mouse)
        continue
    end

    stack = cell_to_string(value_at(rawCells, rowIdx, stackCol));
    address = cell_to_string(value_at(rawCells, rowIdx, addressCol));
    side = normalize_side(cell_to_string(value_at(rawCells, rowIdx, sideCol)));
    inferredSide = infer_side(address);
    if strlength(side) > 0
        lastSide = side;
    elseif strlength(inferredSide) > 0
        side = inferredSide;
        lastSide = side;
    else
        side = lastSide;
    end
    imageFile = filename_from_address(address);

    for regionIdx = 1:numel(cfg.Regions)
        startCol = layerStartCols(regionIdx);
        for markerIdx = 1:cfg.MarkersPerLayer
            marker = cfg.Markers(markerIdx);
            count = cell_to_double(value_at(rawCells, rowIdx, startCol + markerIdx - 1));
            if isnan(count) && cfg.CountRole == "Active"
                count = 0;
            end
            records = append_record(records, cfg, rowIdx, mouse, stack, side, address, ...
                imageFile, cfg.Regions(regionIdx), marker, count);
        end
    end
end
end

function records = import_all_region_totals(rawCells, cfg)
records = empty_record_table();
mouseCol = excel_col_to_num(cfg.MouseCol);
totalCols = arrayfun(@excel_col_to_num, cfg.TotalCols);

for rowIdx = cfg.FirstRow:min(cfg.LastRow, size(rawCells, 1))
    mouse = cell_to_double(value_at(rawCells, rowIdx, mouseCol));
    if isnan(mouse)
        continue
    end

    counts = nan(numel(totalCols), 1);
    for markerIdx = 1:numel(totalCols)
        counts(markerIdx) = cell_to_double(value_at(rawCells, rowIdx, totalCols(markerIdx)));
    end
    if all(isnan(counts))
        continue
    end

    for markerIdx = 1:numel(cfg.Markers)
        count = counts(markerIdx);
        if isnan(count)
            count = nan;
        end
        records = append_record(records, cfg, rowIdx, mouse, "", "", "", "", ...
            "AllRegions", cfg.Markers(markerIdx), count);
    end
end
end

function records = append_record(records, cfg, excelRow, mouse, stack, side, address, ...
    imageFile, region, marker, count)
newRow = table(string(cfg.Dataset), string(cfg.MouseLine), string(cfg.Reporter), ...
    string(cfg.Sheet), string(cfg.MarkerPanel), string(cfg.CountRole), ...
    string(cfg.DenominatorScope), excelRow, mouse, string(stack), string(side), ...
    string(address), string(imageFile), string(region), string(marker), count, ...
    cfg.SlicesRepresented, cfg.FovAreaFactor, ...
    'VariableNames', records.Properties.VariableNames);
records = [records; newRow]; %#ok<AGROW>
end

function records = normalize_record_types(records)
stringVars = {'Dataset', 'MouseLine', 'Reporter', 'Sheet', 'MarkerPanel', ...
    'CountRole', 'DenominatorScope', 'Stack', 'Side', 'Address', ...
    'ImageFile', 'Region', 'Marker'};
for varIdx = 1:numel(stringVars)
    records.(stringVars{varIdx}) = string(records.(stringVars{varIdx}));
end
end

function measurements = make_measurements(records)
activeRows = records.CountRole == "Active";
active = records(activeRows, :);
totals = records(records.CountRole == "Total", :);
inactive = records(records.CountRole == "Inactive", :);

measurements = active(:, {'Dataset', 'MouseLine', 'Reporter', 'Sheet', ...
    'MarkerPanel', 'ExcelRow', 'Mouse', 'Stack', 'Side', 'Address', ...
    'ImageFile', 'Region', 'Marker'});
measurements.Properties.VariableNames{'ExcelRow'} = 'ActiveExcelRow';
measurements.ActiveCount = active.Count;
measurements.ActiveSlicesRepresented = active.SlicesRepresented;
measurements.ActiveFovAreaFactor = active.FovAreaFactor;
measurements.TotalCountRegion = nan(height(active), 1);
measurements.TotalCountAllRegions = nan(height(active), 1);
measurements.DenominatorScope = strings(height(active), 1);
measurements.ActivePercentWithinRegion = nan(height(active), 1);
measurements.ActivePercentOfAllRegions = nan(height(active), 1);

for rowIdx = 1:height(active)
    activeRow = active(rowIdx, :);
    inactiveCount = find_region_inactive_count(inactive, activeRow);
    hasInactiveDenominator = ~isnan(inactiveCount);
    if hasInactiveDenominator
        regionTotal = activeRow.Count + inactiveCount;
    else
        regionTotal = find_region_total(totals, activeRow);
    end

    if hasInactiveDenominator
        allRegionTotal = find_all_regions_from_active_inactive(active, inactive, activeRow);
    else
        allRegionTotal = find_all_region_total(totals, activeRow);
    end

    if isnan(allRegionTotal) && ~isnan(regionTotal)
        allRegionTotal = find_all_regions_from_region_totals(totals, activeRow);
    end

    measurements.TotalCountRegion(rowIdx) = regionTotal;
    measurements.TotalCountAllRegions(rowIdx) = allRegionTotal;
    measurements.ActivePercentWithinRegion(rowIdx) = safe_percent(activeRow.Count, regionTotal);
    measurements.ActivePercentOfAllRegions(rowIdx) = safe_percent(activeRow.Count, allRegionTotal);

    if ~isnan(regionTotal)
        measurements.DenominatorScope(rowIdx) = "region";
    elseif ~isnan(allRegionTotal)
        measurements.DenominatorScope(rowIdx) = "all_regions";
    else
        measurements.DenominatorScope(rowIdx) = "missing";
    end
end
end

function total = find_all_regions_from_active_inactive(active, inactive, activeRow)
activeMask = active.Dataset == activeRow.Dataset & ...
    active.Sheet == activeRow.Sheet & ...
    active.MarkerPanel == activeRow.MarkerPanel & ...
    active.Mouse == activeRow.Mouse & ...
    active.Marker == activeRow.Marker;

activeMask = add_optional_key(activeMask, active.Stack, activeRow.Stack);
activeMask = add_optional_key(activeMask, active.Side, activeRow.Side);
activeMask = add_optional_key(activeMask, active.Address, activeRow.Address);

inactiveMask = inactive.DenominatorScope == "region" & ...
    inactive.Dataset == activeRow.Dataset & ...
    inactive.Sheet == activeRow.Sheet & ...
    inactive.MarkerPanel == activeRow.MarkerPanel & ...
    inactive.Mouse == activeRow.Mouse & ...
    inactive.Marker == activeRow.Marker;

inactiveMask = add_optional_key(inactiveMask, inactive.Stack, activeRow.Stack);
inactiveMask = add_optional_key(inactiveMask, inactive.Side, activeRow.Side);
inactiveMask = add_optional_key(inactiveMask, inactive.Address, activeRow.Address);

activeTotal = sum_or_nan(active.Count(activeMask));
inactiveTotal = sum_or_nan(inactive.Count(inactiveMask));
if isnan(inactiveTotal)
    total = nan;
else
    total = activeTotal + inactiveTotal;
end
end

function count = find_region_inactive_count(inactive, activeRow)
mask = inactive.DenominatorScope == "region" & ...
    inactive.Dataset == activeRow.Dataset & ...
    inactive.Sheet == activeRow.Sheet & ...
    inactive.MarkerPanel == activeRow.MarkerPanel & ...
    inactive.Mouse == activeRow.Mouse & ...
    inactive.Region == activeRow.Region & ...
    inactive.Marker == activeRow.Marker;

mask = add_optional_key(mask, inactive.Stack, activeRow.Stack);
mask = add_optional_key(mask, inactive.Side, activeRow.Side);
count = sum_or_nan(inactive.Count(mask));
end

function total = find_region_total(totals, activeRow)
mask = totals.DenominatorScope == "region" & ...
    totals.Dataset == activeRow.Dataset & ...
    totals.Sheet == activeRow.Sheet & ...
    totals.MarkerPanel == activeRow.MarkerPanel & ...
    totals.Mouse == activeRow.Mouse & ...
    totals.Region == activeRow.Region & ...
    totals.Marker == activeRow.Marker;

mask = add_optional_key(mask, totals.Stack, activeRow.Stack);
mask = add_optional_key(mask, totals.Side, activeRow.Side);
mask = add_optional_key(mask, totals.Address, activeRow.Address);
total = sum_or_nan(totals.Count(mask));
end

function total = find_all_region_total(totals, activeRow)
mask = totals.DenominatorScope == "all_regions" & ...
    totals.Dataset == activeRow.Dataset & ...
    totals.Sheet == activeRow.Sheet & ...
    totals.MarkerPanel == activeRow.MarkerPanel & ...
    totals.Mouse == activeRow.Mouse & ...
    totals.Marker == activeRow.Marker;
total = first_finite_or_nan(totals.Count(mask));
end

function total = find_all_regions_from_region_totals(totals, activeRow)
mask = totals.DenominatorScope == "region" & ...
    totals.Dataset == activeRow.Dataset & ...
    totals.Sheet == activeRow.Sheet & ...
    totals.MarkerPanel == activeRow.MarkerPanel & ...
    totals.Mouse == activeRow.Mouse & ...
    totals.Marker == activeRow.Marker;

mask = add_optional_key(mask, totals.Stack, activeRow.Stack);
mask = add_optional_key(mask, totals.Side, activeRow.Side);
mask = add_optional_key(mask, totals.Address, activeRow.Address);
total = sum_or_nan(totals.Count(mask));
end

function mask = add_optional_key(mask, values, activeValue)
if strlength(activeValue) > 0
    mask = mask & values == activeValue;
end
end

function summary = summarize_measurements(T, groupVars)
summary = empty_summary_table(groupVars);
if isempty(T)
    return
end

hasRegionGrouping = any(strcmp(groupVars, 'Region'));
[G, groupTable] = findgroups(T(:, groupVars));
for groupIdx = 1:height(groupTable)
    rowMask = G == groupIdx;
    groupRows = T(rowMask, :);
    activeValues = groupRows.ActiveCount(isfinite(groupRows.ActiveCount));
    if isempty(activeValues)
        activeCount = nan;
    else
        activeCount = sum(activeValues);
    end

    regionTotals = groupRows.TotalCountRegion(isfinite(groupRows.TotalCountRegion));
    allRegionTotals = groupRows.TotalCountAllRegions(isfinite(groupRows.TotalCountAllRegions));

    if ~isempty(regionTotals)
        totalRegion = sum(regionTotals);
        if hasRegionGrouping && ~isempty(allRegionTotals)
            totalAllRegions = sum(allRegionTotals);
        else
            totalAllRegions = totalRegion;
        end
    else
        totalRegion = nan;
        if hasRegionGrouping
            totalAllRegions = sum_or_nan(allRegionTotals);
        else
            totalAllRegions = first_finite_or_nan(allRegionTotals);
        end
    end

    percentWithinRegion = safe_percent(activeCount, totalRegion);
    percentOfAllRegions = safe_percent(activeCount, totalAllRegions);
    denominatorScope = "missing";
    if ~isnan(totalRegion)
        denominatorScope = "region";
    elseif ~isnan(totalAllRegions)
        denominatorScope = "all_regions";
    end

    nInputRows = height(groupRows);
    nStacks = count_unique_nonempty(groupRows.Stack);
    [nAnalyzedSlices, countNormalizationFactor, meanFovAreaFactor] = ...
        summarize_source_normalization(groupRows);
    nAnalyzedSlices = max(nAnalyzedSlices, 1);
    countNormalizationFactor = max(countNormalizationFactor, 1);
    activeCountPerSlice = activeCount / countNormalizationFactor;
    totalRegionPerSlice = totalRegion / countNormalizationFactor;
    totalAllRegionsPerSlice = totalAllRegions / countNormalizationFactor;
    newRow = [groupTable(groupIdx, :), ...
        table(activeCount, totalRegion, totalAllRegions, percentWithinRegion, ...
            percentOfAllRegions, activeCountPerSlice, totalRegionPerSlice, ...
        totalAllRegionsPerSlice, denominatorScope, nInputRows, nAnalyzedSlices, ...
        meanFovAreaFactor, countNormalizationFactor, nStacks, ...
        'VariableNames', {'ActiveCount', 'TotalCountRegion', ...
        'TotalCountAllRegions', 'ActivePercentWithinRegion', ...
        'ActivePercentOfAllRegions', 'ActiveCountPerSlice', ...
        'TotalCountRegionPerSlice', 'TotalCountAllRegionsPerSlice', ...
        'DenominatorScope', 'NInputRows', 'NAnalyzedSlices', ...
        'MeanFovAreaFactor', 'CountNormalizationFactor', 'NStacks'})];
    summary = [summary; newRow]; %#ok<AGROW>
end
end

function summary = summarize_total_availability(T)
targetRows = ~contains(T.Marker, "-");
T = T(targetRows, :);
if isempty(T)
    summary = table();
    return
end

[G, groupTable] = findgroups(T(:, {'Dataset', 'MouseLine', 'Reporter', 'Sheet', 'MarkerPanel', 'Marker'}));
RowsWithRegionDenominator = nan(height(groupTable), 1);
RowsWithAllRegionDenominator = nan(height(groupTable), 1);
RowsWithMissingDenominator = nan(height(groupTable), 1);
NRows = nan(height(groupTable), 1);

for groupIdx = 1:height(groupTable)
    rows = T(G == groupIdx, :);
    NRows(groupIdx) = height(rows);
    RowsWithRegionDenominator(groupIdx) = sum(rows.DenominatorScope == "region");
    RowsWithAllRegionDenominator(groupIdx) = sum(rows.DenominatorScope == "all_regions");
    RowsWithMissingDenominator(groupIdx) = sum(rows.DenominatorScope == "missing");
end

summary = [groupTable, table(NRows, RowsWithRegionDenominator, ...
    RowsWithAllRegionDenominator, RowsWithMissingDenominator)];
end

function notes = make_notes()
notes = [
    "Datasets are kept separate: cfos_TRE_mKate reports mKate/RFP-tagged cells, intrinsic_cFos reports cFos immunostaining."
    "For intrinsic cFos sheets, region-level denominators are computed as cFos-positive plus cFos-negative marker counts where those reference tables are present."
    "For cfos_TRE_mKate sheets, region-level denominators are computed as RFP/mKate-positive plus RFP/mKate-negative marker counts where those reference tables are present."
    "Mouse-level intrinsic SNCG/CCK8 and proCCK/VIP rows 3:6 and 17:20 use 6 stacks per mouse and a 3x 20x-to-60x field-of-view layer-width normalization."
    "Region-level denominators are imported directly where the workbook contains an irrespective-RFP raw block."
    "Blank cells in total/reference blocks are kept as missing denominators rather than converted to zero."
    "cfos_TRE_mKate slice groups are excluded when they contain no mKate/RFP-positive cells or no positive interneuron-marker cells."
    "S Oriens maps to workbook ABOVE; S Radiatum maps to workbook BELOW."];
end

function data = apply_include_zero_counts(data, includeZeroCounts)
if includeZeroCounts
    return
end
data.records = data.records(data.records.Count > 0, :);
data.measurements = data.measurements(data.measurements.ActiveCount > 0, :);
data.byAnimalRegion = summarize_measurements(data.measurements, ...
    {'Dataset', 'MouseLine', 'Reporter', 'Sheet', 'MarkerPanel', 'Mouse', 'Marker', 'Region'});
data.byAnimalMarker = summarize_measurements(data.measurements, ...
    {'Dataset', 'MouseLine', 'Reporter', 'Sheet', 'MarkerPanel', 'Mouse', 'Marker'});
data.byDatasetMarkerRegion = summarize_measurements(data.measurements, ...
    {'Dataset', 'MouseLine', 'Reporter', 'Sheet', 'MarkerPanel', 'Marker', 'Region'});
data.totalAvailability = summarize_total_availability(data.measurements);
end

function tf = cache_has_required_fields(data)
tf = isfield(data, 'records') && isfield(data, 'measurements') && ...
    isfield(data, 'importerVersion') && data.importerVersion >= 4 && ...
    isfield(data, 'byAnimalRegion') && isfield(data, 'byAnimalMarker') && ...
    isfield(data, 'excludedSlices') && isfield(data, 'excludedRows') && ...
    any(strcmp(data.measurements.Properties.VariableNames, 'ActiveSlicesRepresented')) && ...
    any(strcmp(data.measurements.Properties.VariableNames, 'ActiveFovAreaFactor')) && ...
    any(strcmp(data.byAnimalRegion.Properties.VariableNames, 'ActiveCountPerSlice')) && ...
    any(strcmp(data.byAnimalRegion.Properties.VariableNames, 'NAnalyzedSlices')) && ...
    any(strcmp(data.byAnimalRegion.Properties.VariableNames, 'CountNormalizationFactor'));
end

function data = ensure_exclusion_fields(data)
if ~isfield(data, 'excludedSlices')
    data.excludedSlices = empty_excluded_slice_table();
end
if ~isfield(data, 'excludedRows')
    if isfield(data, 'records') && istable(data.records)
        data.excludedRows = data.records(false(height(data.records), 1), :);
        data.excludedRows.ExclusionReason = strings(0, 1);
    else
        data.excludedRows = table();
    end
end
end

function [records, excludedSlices, excludedRows] = exclude_empty_mkate_slices(records)
excludedSlices = empty_excluded_slice_table();
excludedRows = records(false(height(records), 1), :);
excludedRows.ExclusionReason = strings(0, 1);

if isempty(records)
    return
end

sliceVars = {'Dataset', 'MouseLine', 'Reporter', 'Sheet', 'MarkerPanel', ...
    'Mouse', 'Stack', 'Side', 'Address', 'ImageFile'};
[G, groupTable] = findgroups(records(:, sliceVars));
excludeMask = false(height(records), 1);
exclusionReasons = strings(height(records), 1);

for groupIdx = 1:height(groupTable)
    rowMask = G == groupIdx;
    rows = records(rowMask, :);
    if rows.Dataset(1) ~= "cfos_TRE_mKate" || ~any(rows.CountRole == "Active")
        continue
    end

    mKatePositiveCount = sum_positive_counts(rows.Count(rows.CountRole == "Active"));
    markerPositiveCount = marker_positive_total_count(rows);
    reasons = strings(0, 1);
    if mKatePositiveCount <= 0
        reasons(end + 1, 1) = "No mKate/RFP-positive cells in slice"; %#ok<AGROW>
    end
    if markerPositiveCount <= 0
        reasons(end + 1, 1) = "No positive interneuron-marker cells in slice"; %#ok<AGROW>
    end

    if isempty(reasons)
        continue
    end

    reason = strjoin(reasons, "; ");
    excludeMask(rowMask) = true;
    exclusionReasons(rowMask) = reason;
    excludedSlices = [excludedSlices; make_excluded_slice_row( ... %#ok<AGROW>
        groupTable(groupIdx, :), rows, mKatePositiveCount, markerPositiveCount, reason)];
end

excludedRows = records(excludeMask, :);
excludedRows.ExclusionReason = exclusionReasons(excludeMask);
records = records(~excludeMask, :);
end

function count = marker_positive_total_count(rows)
positiveMarkerRows = is_positive_marker(rows.Marker);
totalRows = rows.CountRole == "Total" & positiveMarkerRows;
if any(isfinite(rows.Count(totalRows)))
    count = sum_positive_counts(rows.Count(totalRows));
    return
end

referenceRows = (rows.CountRole == "Active" | rows.CountRole == "Inactive") & ...
    positiveMarkerRows;
count = sum_positive_counts(rows.Count(referenceRows));
end

function tf = is_positive_marker(markers)
tf = strlength(markers) > 0 & ~contains(markers, "-");
end

function count = sum_positive_counts(values)
values = values(isfinite(values) & values > 0);
if isempty(values)
    count = 0;
else
    count = sum(values);
end
end

function row = make_excluded_slice_row(groupValues, rows, mKatePositiveCount, ...
    markerPositiveCount, reason)
row = table(groupValues.Dataset, groupValues.MouseLine, groupValues.Reporter, ...
    groupValues.Sheet, groupValues.MarkerPanel, groupValues.Mouse, ...
    groupValues.Stack, groupValues.Side, groupValues.Address, groupValues.ImageFile, ...
    mKatePositiveCount, markerPositiveCount, height(rows), ...
    format_excel_rows(rows.ExcelRow(rows.CountRole == "Active")), ...
    format_excel_rows(rows.ExcelRow(rows.CountRole == "Total")), ...
    format_excel_rows(rows.ExcelRow(rows.CountRole == "Inactive")), ...
    string(reason), 'VariableNames', empty_excluded_slice_variable_names());
end

function label = format_excel_rows(rows)
rows = unique(rows(isfinite(rows)));
if isempty(rows)
    label = "";
else
    label = strjoin(string(rows(:)'), ", ");
end
end

function T = empty_excluded_slice_table()
T = table(strings(0, 1), strings(0, 1), strings(0, 1), strings(0, 1), ...
    strings(0, 1), nan(0, 1), strings(0, 1), strings(0, 1), strings(0, 1), ...
    strings(0, 1), nan(0, 1), nan(0, 1), nan(0, 1), strings(0, 1), ...
    strings(0, 1), strings(0, 1), strings(0, 1), ...
    'VariableNames', empty_excluded_slice_variable_names());
end

function names = empty_excluded_slice_variable_names()
names = {'Dataset', 'MouseLine', 'Reporter', 'Sheet', 'MarkerPanel', ...
    'Mouse', 'Stack', 'Side', 'Address', 'ImageFile', 'MKatePositiveCount', ...
    'InterneuronMarkerPositiveCount', 'ExcludedRecordCount', ...
    'ActiveExcelRows', 'TotalExcelRows', 'InactiveExcelRows', ...
    'ExclusionReason'};
end

function records = empty_record_table()
records = cell2table(cell(0, 18), 'VariableNames', ...
    {'Dataset', 'MouseLine', 'Reporter', 'Sheet', 'MarkerPanel', ...
    'CountRole', 'DenominatorScope', 'ExcelRow', 'Mouse', 'Stack', ...
    'Side', 'Address', 'ImageFile', 'Region', 'Marker', 'Count', ...
    'SlicesRepresented', 'FovAreaFactor'});
records.Dataset = strings(0, 1);
records.MouseLine = strings(0, 1);
records.Reporter = strings(0, 1);
records.Sheet = strings(0, 1);
records.MarkerPanel = strings(0, 1);
records.CountRole = strings(0, 1);
records.DenominatorScope = strings(0, 1);
records.ExcelRow = nan(0, 1);
records.Mouse = nan(0, 1);
records.Stack = strings(0, 1);
records.Side = strings(0, 1);
records.Address = strings(0, 1);
records.ImageFile = strings(0, 1);
records.Region = strings(0, 1);
records.Marker = strings(0, 1);
records.Count = nan(0, 1);
records.SlicesRepresented = nan(0, 1);
records.FovAreaFactor = nan(0, 1);
end

function summary = empty_summary_table(groupVars)
summary = table();
for varIdx = 1:numel(groupVars)
    varName = groupVars{varIdx};
    if strcmp(varName, 'Mouse')
        summary.(varName) = nan(0, 1);
    else
        summary.(varName) = strings(0, 1);
    end
end
summary.ActiveCount = nan(0, 1);
summary.TotalCountRegion = nan(0, 1);
summary.TotalCountAllRegions = nan(0, 1);
summary.ActivePercentWithinRegion = nan(0, 1);
summary.ActivePercentOfAllRegions = nan(0, 1);
summary.ActiveCountPerSlice = nan(0, 1);
summary.TotalCountRegionPerSlice = nan(0, 1);
summary.TotalCountAllRegionsPerSlice = nan(0, 1);
summary.DenominatorScope = strings(0, 1);
summary.NInputRows = nan(0, 1);
summary.NAnalyzedSlices = nan(0, 1);
summary.MeanFovAreaFactor = nan(0, 1);
summary.CountNormalizationFactor = nan(0, 1);
summary.NStacks = nan(0, 1);
end

function value = value_at(cells, rowIdx, colIdx)
if isnan(colIdx) || rowIdx < 1 || colIdx < 1 || rowIdx > size(cells, 1) || colIdx > size(cells, 2)
    value = [];
else
    value = cells{rowIdx, colIdx};
end
end

function value = cell_to_double(cellValue)
value = nan;
if isempty(cellValue)
    return
elseif isnumeric(cellValue) && isscalar(cellValue)
    value = double(cellValue);
elseif islogical(cellValue) && isscalar(cellValue)
    value = double(cellValue);
elseif ischar(cellValue) || isstring(cellValue)
    textValue = strtrim(string(cellValue));
    if strlength(textValue) > 0 && ~ismissing(textValue)
        parsed = str2double(textValue);
        if ~isnan(parsed)
            value = parsed;
        end
    end
end
end

function value = cell_to_string(cellValue)
if isempty(cellValue)
    value = "";
elseif isstring(cellValue)
    if ismissing(cellValue)
        value = "";
    else
        value = string(cellValue);
    end
elseif ischar(cellValue)
    value = string(cellValue);
elseif isnumeric(cellValue) || islogical(cellValue)
    if isscalar(cellValue) && isnan(cellValue)
        value = "";
    else
        value = string(cellValue);
    end
else
    value = string(cellValue);
end
value = strtrim(value);
end

function side = normalize_side(sideText)
textValue = lower(strtrim(string(sideText)));
if strlength(textValue) == 0 || ismissing(textValue)
    side = "";
elseif contains(textValue, "non")
    side = "Non-implanted";
elseif contains(textValue, "impl")
    side = "Implanted";
else
    side = strtrim(string(sideText));
end
end

function side = infer_side(address)
textValue = lower(string(address));
if strlength(textValue) == 0
    side = "";
elseif contains(textValue, "non-implanted") || contains(textValue, "non_implanted") || ...
        contains(textValue, "nonimplant") || contains(textValue, "non-implant") || ...
        contains(textValue, "nonimplant")
    side = "Non-implanted";
elseif contains(textValue, "implanted") || contains(textValue, "implant") || ...
        contains(textValue, "inj")
    side = "Implanted";
else
    side = "";
end
end

function imageFile = filename_from_address(address)
if strlength(address) == 0
    imageFile = "";
    return
end
parts = regexp(char(address), '[\\/]', 'split');
imageFile = string(parts{end});
end

function idx = excel_col_to_num(col)
if strlength(string(col)) == 0
    idx = nan;
    return
end
col = char(upper(string(col)));
idx = 0;
for charIdx = 1:numel(col)
    idx = idx * 26 + double(col(charIdx)) - double('A') + 1;
end
end

function value = safe_percent(numerator, denominator)
if isnan(denominator) || denominator == 0
    value = nan;
else
    value = numerator * 100 / denominator;
end
end

function value = sum_or_nan(values)
values = values(isfinite(values));
if isempty(values)
    value = nan;
else
    value = sum(values);
end
end

function value = first_finite_or_nan(values)
values = values(isfinite(values));
if isempty(values)
    value = nan;
else
    value = values(1);
end
end

function n = count_unique_nonempty(values)
values = values(strlength(values) > 0);
if isempty(values)
    n = 0;
else
    n = numel(unique(values));
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

function [nSlices, normalizationFactor, meanFovAreaFactor] = summarize_source_normalization(rows)
nSlices = 0;
normalizationFactor = 0;
meanFovAreaFactor = 1;

if isempty(rows)
    return
end

hasSlices = any(strcmp(rows.Properties.VariableNames, 'ActiveSlicesRepresented'));
hasFovFactor = any(strcmp(rows.Properties.VariableNames, 'ActiveFovAreaFactor'));
if ~hasSlices
    nSlices = count_unique_finite(rows.ActiveExcelRow);
    normalizationFactor = nSlices;
    return
end

sourceRows = rows(isfinite(rows.ActiveExcelRow), :);
if isempty(sourceRows)
    return
end

[~, uniqueIdx] = unique(sourceRows.ActiveExcelRow, 'stable');
sliceValues = sourceRows.ActiveSlicesRepresented(uniqueIdx);
sliceValues(~isfinite(sliceValues) | sliceValues <= 0) = 1;

if hasFovFactor
    fovValues = sourceRows.ActiveFovAreaFactor(uniqueIdx);
else
    fovValues = ones(size(sliceValues));
end
fovValues(~isfinite(fovValues) | fovValues <= 0) = 1;

nSlices = sum(sliceValues);
normalizationFactor = sum(sliceValues .* fovValues);
if nSlices > 0
    meanFovAreaFactor = normalizationFactor / nSlices;
end
end

function workbookFile = default_workbook_file()
codeFolder = fileparts(mfilename('fullpath'));
repoRoot = fullfile(codeFolder, '..', '..');
candidateFiles = {
    fullfile(repoRoot, 'Results', 'Quantifications.xlsx')
    fullfile(repoRoot, 'SDC microscopy', 'Quantifications.xlsx')
    fullfile(repoRoot, '..', 'Results', 'Quantifications.xlsx')
    fullfile(repoRoot, '..', 'SDC microscopy', 'Quantifications.xlsx')};
workbookFile = first_existing_file(candidateFiles);
end

function cacheFile = default_cache_file(workbookFile)
cacheFile = fullfile(fileparts(char(workbookFile)), 'interneuron_quantifications_import.mat');
end

function filePath = first_existing_file(candidateFiles)
filePath = candidateFiles{1};
for idx = 1:numel(candidateFiles)
    if isfile(candidateFiles{idx})
        filePath = candidateFiles{idx};
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
validNames = {'includezerocounts', 'cachefile', 'reimportexcel'};
tf = any(strcmp(name, validNames));
end
