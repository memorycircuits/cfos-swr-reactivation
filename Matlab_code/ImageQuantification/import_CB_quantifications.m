function data = import_CB_quantifications(workbookFile, varargin)
%IMPORT_CB_QUANTIFICATIONS Load CB/RFP quantifications, preferring the MAT cache.

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
addParameter(parser, 'Sheet', 'CB+RFP', @(x) ischar(x) || isstring(x));
addParameter(parser, 'MKateSheet', 'MKate', @(x) ischar(x) || isstring(x));
addParameter(parser, 'IncludeZeroCounts', true, @(x) islogical(x) && isscalar(x));
addParameter(parser, 'CacheFile', '', @(x) ischar(x) || isstring(x));
addParameter(parser, 'ReimportExcel', false, @(x) islogical(x) && isscalar(x));
parse(parser, varargin{:});

sheetName = char(parser.Results.Sheet);
mkateSheetName = char(parser.Results.MKateSheet);
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
        data.loadedFrom = "mat";
        data = apply_include_zero_counts(data, includeZeroCounts);
        return
    end
end

data = import_cb_sheet_from_excel(workbookFile, sheetName, mkateSheetName);
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

function data = import_cb_sheet_from_excel(workbookFile, sheetName, mkateSheetName)
if ~isfile(workbookFile)
    error('CBImport:MissingWorkbook', 'Workbook not found: %s', char(workbookFile));
end

availableSheets = sheetnames(workbookFile);
if ~any(strcmp(availableSheets, sheetName))
    error('CBImport:MissingSheet', 'Sheet "%s" not found in %s.', sheetName, char(workbookFile));
end

rawCells = readcell(workbookFile, 'Sheet', sheetName, 'Range', 'A:AQ');
rowNumbers = (1:size(rawCells, 1))';
addresses = cells_to_string(rawCells(:, 3));
dataRowMask = rowNumbers > 3 & strlength(strtrim(addresses)) > 0;
dataRows = rowNumbers(dataRowMask);

countColumns = 4:15;
counts = zeros(numel(dataRows), numel(countColumns));
for colIdx = 1:numel(countColumns)
    colValues = cells_to_double(rawCells(dataRowMask, countColumns(colIdx)));
    colValues(isnan(colValues)) = 0;
    counts(:, colIdx) = colValues;
end

ExcelRow = dataRows;
Mouse = cells_to_double(rawCells(dataRowMask, 1));
Stack = cells_to_string(rawCells(dataRowMask, 2));
Address = addresses(dataRowMask);
Side = strings(numel(dataRows), 1);
ImageFile = strings(numel(dataRows), 1);
for rowIdx = 1:numel(dataRows)
    Side(rowIdx) = infer_side(Address(rowIdx));
    ImageFile(rowIdx) = filename_from_address(Address(rowIdx));
end

rawTable = table(ExcelRow, Mouse, Stack, Side, Address, ImageFile);
countNames = [
    "RfpCbPositiveSuperficial"
    "RfpCbPositiveDeep"
    "RfpCbPositiveOriens"
    "RfpCbPositiveRadiatum"
    "RfpCbNegativeSuperficial"
    "RfpCbNegativeDeep"
    "RfpCbNegativeOriens"
    "RfpCbNegativeRadiatum"
    "RfpCbAmbiguousSuperficial"
    "RfpCbAmbiguousDeep"
    "RfpCbAmbiguousOriens"
    "RfpCbAmbiguousRadiatum"];

for colIdx = 1:numel(countNames)
    rawTable.(char(countNames(colIdx))) = counts(:, colIdx);
end

rawTable.TotalRfp = sum(counts, 2);
rawTable.TotalSuperficial = sum(counts(:, [1 5 9]), 2);
rawTable.TotalDeep = sum(counts(:, [2 6 10]), 2);
rawTable.TotalOriens = sum(counts(:, [3 7 11]), 2);
rawTable.TotalRadiatum = sum(counts(:, [4 8 12]), 2);
rawTable.TotalCbPositive = sum(counts(:, 1:4), 2);
rawTable.TotalCbNegative = sum(counts(:, 5:8), 2);
rawTable.TotalCbAmbiguous = sum(counts(:, 9:12), 2);
rawTable.CbPositiveSuperficial = counts(:, 1);
rawTable.CbPositiveDeep = counts(:, 2);
rawTable.CbPositiveOriens = counts(:, 3);
rawTable.CbPositiveRadiatum = counts(:, 4);
rawTable.CbNegativeSuperficial = counts(:, 5);
rawTable.CbNegativeDeep = counts(:, 6);
rawTable.CbNegativeOriens = counts(:, 7);
rawTable.CbNegativeRadiatum = counts(:, 8);
rawTable.CbAmbiguousSuperficial = counts(:, 9);
rawTable.CbAmbiguousDeep = counts(:, 10);
rawTable.CbAmbiguousOriens = counts(:, 11);
rawTable.CbAmbiguousRadiatum = counts(:, 12);

longTable = make_long_table(rawTable, counts, true);

data = struct();
data.workbookFile = string(workbookFile);
data.sheetName = string(sheetName);
data.raw = rawTable;
data.long = longTable;
data.mKate = import_mkate_sheet_from_excel(workbookFile, availableSheets, mkateSheetName);
data = add_summary_tables(data);

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
cacheFile = fullfile(fileparts(char(workbookFile)), 'CB_quantifications_import.mat');
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

function tf = cache_has_required_fields(data)
tf = isfield(data, 'raw') && isfield(data, 'long') && ...
    isfield(data, 'mKate') && isfield(data.mKate, 'layerAllocation') && ...
    ~isempty(data.mKate.layerAllocation) && height(data.mKate.layerAllocation) > 0;
end

function mKate = import_mkate_sheet_from_excel(workbookFile, availableSheets, requestedSheetName)
mKate = struct();
mKate.requestedSheetName = string(requestedSheetName);
mKate.sheetName = "";
mKate.raw = table();
mKate.layerAllocation = table();

resolvedSheetName = resolve_sheet_name(availableSheets, requestedSheetName);
if strlength(resolvedSheetName) == 0
    warning('CBImport:MissingMKateSheet', ...
        'mKate sheet "%s" not found in %s. Panel 1 will use the CB+RFP sheet.', ...
        requestedSheetName, char(workbookFile));
    return
end

rawCells = readcell(workbookFile, 'Sheet', char(resolvedSheetName), 'Range', 'A:O');
rowNumbers = (1:size(rawCells, 1))';
mouseValues = cells_to_double(rawCells(:, 1));
dataRowMask = rowNumbers > 1 & ~isnan(mouseValues);
dataRows = rowNumbers(dataRowMask);

if isempty(dataRows)
    warning('CBImport:EmptyMKateSheet', ...
        'mKate sheet "%s" did not contain numeric mouse rows. Panel 1 will use the CB+RFP sheet.', ...
        char(resolvedSheetName));
    return
end

countColumns = 2:5;
counts = zeros(numel(dataRows), numel(countColumns));
for colIdx = 1:numel(countColumns)
    colValues = cells_to_double(rawCells(dataRowMask, countColumns(colIdx)));
    colValues(isnan(colValues)) = 0;
    counts(:, colIdx) = colValues;
end

ExcelRow = dataRows;
Mouse = mouseValues(dataRowMask);
AboveCount = counts(:, 1);
DeepCount = counts(:, 2);
SuperficialCount = counts(:, 3);
BelowCount = counts(:, 4);
TotalMKateCount = sum(counts, 2);
InLayerTotal = SuperficialCount + DeepCount;
SuperficialPercentOfInLayer = nan(numel(dataRows), 1);
DeepPercentOfInLayer = nan(numel(dataRows), 1);
DeepMinusSuperficialPercent = nan(numel(dataRows), 1);

for rowIdx = 1:numel(dataRows)
    SuperficialPercentOfInLayer(rowIdx) = safe_percent(SuperficialCount(rowIdx), InLayerTotal(rowIdx));
    DeepPercentOfInLayer(rowIdx) = safe_percent(DeepCount(rowIdx), InLayerTotal(rowIdx));
    DeepMinusSuperficialPercent(rowIdx) = DeepPercentOfInLayer(rowIdx) - SuperficialPercentOfInLayer(rowIdx);
end

rawTable = table(ExcelRow, Mouse, AboveCount, DeepCount, SuperficialCount, ...
    BelowCount, TotalMKateCount, InLayerTotal, SuperficialPercentOfInLayer, ...
    DeepPercentOfInLayer, DeepMinusSuperficialPercent);

mKate.sheetName = resolvedSheetName;
mKate.raw = rawTable;
mKate.layerAllocation = rawTable(:, {'Mouse', 'SuperficialCount', 'DeepCount', ...
    'InLayerTotal', 'SuperficialPercentOfInLayer', 'DeepPercentOfInLayer', ...
    'DeepMinusSuperficialPercent', 'AboveCount', 'BelowCount', 'TotalMKateCount'});
end

function sheetName = resolve_sheet_name(availableSheets, requestedSheetName)
requested = string(requestedSheetName);
available = string(availableSheets);
matchIdx = find(strcmpi(available, requested), 1);
if isempty(matchIdx)
    sheetName = "";
else
    sheetName = available(matchIdx);
end
end

function data = apply_include_zero_counts(data, includeZeroCounts)
if includeZeroCounts
    return
end
data.long = data.long(data.long.Count > 0, :);
data = add_summary_tables(data);
end

function data = add_summary_tables(data)
longTable = data.long;
data.byMouseStatusPosition = rename_summary_count(groupsummary(longTable, ...
    {'Mouse', 'CbStatus', 'Position'}, 'sum', 'Count'));
data.byMouseStatus = rename_summary_count(groupsummary(longTable, ...
    {'Mouse', 'CbStatus'}, 'sum', 'Count'));
data.byMousePosition = rename_summary_count(groupsummary(longTable, ...
    {'Mouse', 'Position'}, 'sum', 'Count'));
data.byMouseSideStatusPosition = rename_summary_count(groupsummary(longTable, ...
    {'Mouse', 'Side', 'CbStatus', 'Position'}, 'sum', 'Count'));
end

function longTable = make_long_table(rawTable, counts, includeZeroCounts)
statusNames = ["CB+"; "CB-"; "CB+/-"];
positionNames = ["Superficial"; "Deep"; "Oriens"; "Radiatum"];
layerBins = ["InLayer"; "InLayer"; "OutsideLayer"; "OutsideLayer"];
statusOffsets = [0 4 8];

nRows = height(rawTable);
nLong = nRows * numel(statusNames) * numel(positionNames);
ExcelRow = nan(nLong, 1);
Mouse = nan(nLong, 1);
Stack = strings(nLong, 1);
Side = strings(nLong, 1);
Address = strings(nLong, 1);
ImageFile = strings(nLong, 1);
CbStatus = strings(nLong, 1);
Position = strings(nLong, 1);
LayerBin = strings(nLong, 1);
Count = nan(nLong, 1);

k = 0;
for rowIdx = 1:nRows
    for statusIdx = 1:numel(statusNames)
        for positionIdx = 1:numel(positionNames)
            k = k + 1;
            ExcelRow(k) = rawTable.ExcelRow(rowIdx);
            Mouse(k) = rawTable.Mouse(rowIdx);
            Stack(k) = rawTable.Stack(rowIdx);
            Side(k) = rawTable.Side(rowIdx);
            Address(k) = rawTable.Address(rowIdx);
            ImageFile(k) = rawTable.ImageFile(rowIdx);
            CbStatus(k) = statusNames(statusIdx);
            Position(k) = positionNames(positionIdx);
            LayerBin(k) = layerBins(positionIdx);
            Count(k) = counts(rowIdx, statusOffsets(statusIdx) + positionIdx);
        end
    end
end

longTable = table(ExcelRow, Mouse, Stack, Side, Address, ImageFile, ...
    CbStatus, Position, LayerBin, Count);

if ~includeZeroCounts
    longTable = longTable(longTable.Count > 0, :);
end
end

function outTable = rename_summary_count(inTable)
outTable = inTable;
names = outTable.Properties.VariableNames;
idx = strcmp(names, 'sum_Count');
if any(idx)
    names{idx} = 'Count';
    outTable.Properties.VariableNames = names;
end
end

function values = cells_to_string(cells)
values = strings(numel(cells), 1);
for idx = 1:numel(cells)
    value = cells{idx};
    if is_empty_cell(value)
        values(idx) = "";
    elseif isstring(value)
        if ismissing(value)
            values(idx) = "";
        else
            values(idx) = value;
        end
    elseif ischar(value)
        values(idx) = string(value);
    elseif isnumeric(value) || islogical(value)
        values(idx) = string(value);
    else
        values(idx) = string(value);
    end
end
end

function values = cells_to_double(cells)
values = nan(numel(cells), 1);
for idx = 1:numel(cells)
    value = cells{idx};
    if is_empty_cell(value)
        continue
    elseif isnumeric(value) && isscalar(value)
        values(idx) = double(value);
    elseif islogical(value) && isscalar(value)
        values(idx) = double(value);
    elseif ischar(value) || isstring(value)
        textValue = strtrim(string(value));
        if strlength(textValue) > 0 && ~ismissing(textValue)
            values(idx) = str2double(textValue);
        end
    end
end
end

function value = safe_percent(numerator, denominator)
if denominator == 0
    value = nan;
else
    value = numerator * 100 / denominator;
end
end

function tf = is_empty_cell(value)
tf = isempty(value) || ...
    (isnumeric(value) && isscalar(value) && isnan(value)) || ...
    (isstring(value) && (isempty(value) || ismissing(value) || strlength(value) == 0)) || ...
    (ischar(value) && isempty(strtrim(value)));
end

function side = infer_side(address)
textValue = lower(string(address));
if contains(textValue, "non-implanted") || contains(textValue, "non_inj") || ...
        contains(textValue, "non-inj") || contains(textValue, "non-injection") || ...
        contains(textValue, "non injection")
    side = "Non-implanted";
elseif contains(textValue, "implanted side") || contains(textValue, "injection site") || ...
        contains(textValue, "_inj") || contains(textValue, "-inj")
    side = "Implanted";
else
    side = "Unknown";
end
end

function imageFile = filename_from_address(address)
parts = regexp(char(address), '[\\/]', 'split');
imageFile = string(parts{end});
end

function tf = is_name_value_start(value)
if ~(ischar(value) || (isstring(value) && isscalar(value)))
    tf = false;
    return
end
name = lower(char(string(value)));
validNames = {'sheet', 'mkatesheet', 'includezerocounts', 'cachefile', 'reimportexcel'};
tf = any(strcmp(name, validNames));
end
