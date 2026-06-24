function report = check_interneuron_quantification_math(workbookFile, varargin)
%CHECK_INTERNEURON_QUANTIFICATION_MATH Check interneuron workbook consistency.
%   REPORT = CHECK_INTERNEURON_QUANTIFICATION_MATH(WORKBOOKFILE) reimports
%   Quantifications.xlsx and checks for formula errors, invalid count
%   values, active counts exceeding denominators, and missing denominators.
%
%   Name-value inputs:
%     Tolerance  Numeric tolerance for integer and denominator checks.
%     CacheFile  MAT cache path passed to the importer.
%
%   REPORT is a table with one row per check and a Passed flag. Calling the
%   function without an output prints failed checks and denominator coverage.

if nargin < 1
    workbookFile = [];
elseif is_name_value_start(workbookFile)
    varargin = [{workbookFile}, varargin];
    workbookFile = [];
end

parser = inputParser;
addParameter(parser, 'Tolerance', 1e-9, @(x) isnumeric(x) && isscalar(x) && x >= 0);
addParameter(parser, 'CacheFile', '', @(x) ischar(x) || isstring(x));
parse(parser, varargin{:});

tolerance = parser.Results.Tolerance;
cacheFile = parser.Results.CacheFile;

data = import_interneuron_quantifications(workbookFile, 'IncludeZeroCounts', true, ...
    'CacheFile', cacheFile, 'ReimportExcel', true);

categories = strings(0, 1);
details = strings(0, 1);
datasets = strings(0, 1);
sheets = strings(0, 1);
markerPanels = strings(0, 1);
mice = nan(0, 1);
sides = strings(0, 1);
stacks = strings(0, 1);
regions = strings(0, 1);
markers = strings(0, 1);
activeCounts = nan(0, 1);
denominators = nan(0, 1);
differences = nan(0, 1);
passed = false(0, 1);

check_formula_errors(data.workbookFile);
check_count_values(data.records);
check_region_denominators(data.measurements);
check_summary_region_denominators(data.byAnimalRegion);
check_all_region_denominators(data.byAnimalMarker);
check_missing_denominators(data.measurements);

report = table(categories, details, datasets, sheets, markerPanels, mice, sides, ...
    stacks, regions, markers, activeCounts, denominators, differences, passed, ...
    'VariableNames', {'Category', 'Detail', 'Dataset', 'Sheet', 'MarkerPanel', ...
    'Mouse', 'Side', 'Stack', 'Region', 'Marker', 'ActiveCount', ...
    'Denominator', 'Difference', 'Passed'});

if nargout == 0
    failed = report(~report.Passed, :);
    if isempty(failed)
        fprintf('All %d interneuron workbook checks passed within tolerance %.3g.\n', ...
            height(report), tolerance);
    else
        fprintf('%d of %d interneuron workbook checks failed within tolerance %.3g.\n', ...
            height(failed), height(report), tolerance);
        disp(failed(:, {'Category', 'Detail', 'Dataset', 'Sheet', 'MarkerPanel', ...
            'Mouse', 'Side', 'Stack', 'Region', 'Marker', 'ActiveCount', ...
            'Denominator', 'Difference'}));
    end

    if isfield(data, 'totalAvailability') && ~isempty(data.totalAvailability)
        fprintf('\nDenominator availability by dataset/marker:\n');
        disp(data.totalAvailability);
    end
end

    function check_formula_errors(workbookPath)
        scanSheets = ["PV+SOM+RFP", "PV+CFOS+SOM", "aSNCG+cFos+CCK8"];
        availableSheets = string(sheetnames(workbookPath));
        for sheetIdx = 1:numel(scanSheets)
            sheetName = scanSheets(sheetIdx);
            if ~any(strcmpi(availableSheets, sheetName))
                continue
            end
            rawCells = readcell(workbookPath, 'Sheet', char(sheetName), 'Range', 'A:DX');
            for rowIdx = 1:size(rawCells, 1)
                for colIdx = 1:size(rawCells, 2)
                    value = rawCells{rowIdx, colIdx};
                    if is_formula_error_value(value)
                        add_check('formula error', sprintf('%s%d contains %s', ...
                            char(excel_num_to_col(colIdx)), rowIdx, char(string(value))), ...
                            "", sheetName, "", nan, "", "", "", "", nan, nan, nan, false);
                    end
                end
            end
        end
    end

    function check_count_values(records)
        for rowIdx = 1:height(records)
            count = records.Count(rowIdx);
            if isnan(count)
                continue
            end
            isValid = count >= -tolerance && abs(count - round(count)) <= tolerance;
            add_check('count value', 'Counts should be non-negative integers.', ...
                records.Dataset(rowIdx), records.Sheet(rowIdx), records.MarkerPanel(rowIdx), ...
                records.Mouse(rowIdx), records.Side(rowIdx), records.Stack(rowIdx), ...
                records.Region(rowIdx), records.Marker(rowIdx), count, nan, nan, isValid);
        end
    end

    function check_region_denominators(measurements)
        rows = isfinite(measurements.TotalCountRegion);
        for rowIdx = find(rows)'
            active = measurements.ActiveCount(rowIdx);
            denom = measurements.TotalCountRegion(rowIdx);
            diffValue = active - denom;
            isValid = active <= denom + tolerance;
            add_check('active <= region total', ...
                'Active marker count should not exceed the region-level marker total.', ...
                measurements.Dataset(rowIdx), measurements.Sheet(rowIdx), measurements.MarkerPanel(rowIdx), ...
                measurements.Mouse(rowIdx), measurements.Side(rowIdx), measurements.Stack(rowIdx), ...
                measurements.Region(rowIdx), measurements.Marker(rowIdx), active, denom, diffValue, isValid);
        end
    end

    function check_all_region_denominators(summaryTable)
        rows = isfinite(summaryTable.TotalCountAllRegions);
        for rowIdx = find(rows)'
            active = summaryTable.ActiveCount(rowIdx);
            denom = summaryTable.TotalCountAllRegions(rowIdx);
            diffValue = active - denom;
            isValid = active <= denom + tolerance;
            add_check('active <= all-region total', ...
                'Animal-level active count should not exceed the all-region marker total.', ...
                summaryTable.Dataset(rowIdx), summaryTable.Sheet(rowIdx), summaryTable.MarkerPanel(rowIdx), ...
                summaryTable.Mouse(rowIdx), "", "", "", summaryTable.Marker(rowIdx), ...
                active, denom, diffValue, isValid);
        end
    end

    function check_summary_region_denominators(summaryTable)
        rows = isfinite(summaryTable.TotalCountRegion);
        for rowIdx = find(rows)'
            active = summaryTable.ActiveCount(rowIdx);
            denom = summaryTable.TotalCountRegion(rowIdx);
            diffValue = active - denom;
            isValid = active <= denom + tolerance;
            add_check('summary active <= region total', ...
                'Animal-region active count should not exceed the summarized region-level marker total.', ...
                summaryTable.Dataset(rowIdx), summaryTable.Sheet(rowIdx), summaryTable.MarkerPanel(rowIdx), ...
                summaryTable.Mouse(rowIdx), "", "", summaryTable.Region(rowIdx), ...
                summaryTable.Marker(rowIdx), active, denom, diffValue, isValid);
        end
    end

    function check_missing_denominators(measurements)
        targetRows = measurements.DenominatorScope == "missing" & ...
            ~contains(measurements.Marker, "-");
        for rowIdx = find(targetRows)'
            add_check('denominator availability', ...
                'No matching total marker count was imported for this active target-marker row.', ...
                measurements.Dataset(rowIdx), measurements.Sheet(rowIdx), measurements.MarkerPanel(rowIdx), ...
                measurements.Mouse(rowIdx), measurements.Side(rowIdx), measurements.Stack(rowIdx), ...
                measurements.Region(rowIdx), measurements.Marker(rowIdx), ...
                measurements.ActiveCount(rowIdx), nan, nan, false);
        end
    end

    function add_check(category, detail, dataset, sheet, markerPanel, mouse, side, ...
            stack, region, marker, activeCount, denominator, difference, didPass)
        categories(end + 1, 1) = string(category);
        details(end + 1, 1) = string(detail);
        datasets(end + 1, 1) = string(dataset);
        sheets(end + 1, 1) = string(sheet);
        markerPanels(end + 1, 1) = string(markerPanel);
        mice(end + 1, 1) = mouse;
        sides(end + 1, 1) = string(side);
        stacks(end + 1, 1) = string(stack);
        regions(end + 1, 1) = string(region);
        markers(end + 1, 1) = string(marker);
        activeCounts(end + 1, 1) = activeCount;
        denominators(end + 1, 1) = denominator;
        differences(end + 1, 1) = difference;
        passed(end + 1, 1) = didPass;
    end
end

function tf = is_formula_error_value(value)
tf = false;
if ischar(value) || isstring(value)
    textValue = strtrim(string(value));
    tf = startsWith(textValue, "#");
end
end

function col = excel_num_to_col(idx)
chars = '';
while idx > 0
    remValue = mod(idx - 1, 26);
    chars = [char(double('A') + remValue), chars]; %#ok<AGROW>
    idx = floor((idx - 1) / 26);
end
col = string(chars);
end

function tf = is_name_value_start(value)
if ~(ischar(value) || (isstring(value) && isscalar(value)))
    tf = false;
    return
end
name = lower(char(string(value)));
validNames = {'tolerance', 'cachefile'};
tf = any(strcmp(name, validNames));
end
