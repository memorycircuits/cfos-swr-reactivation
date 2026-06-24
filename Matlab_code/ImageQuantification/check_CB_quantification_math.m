function report = check_CB_quantification_math(workbookFile, varargin)
%CHECK_CB_QUANTIFICATION_MATH Recompute CB worksheet totals from raw counts.

if nargin < 1
    workbookFile = [];
elseif is_name_value_start(workbookFile)
    varargin = [{workbookFile}, varargin];
    workbookFile = [];
end

parser = inputParser;
addParameter(parser, 'Sheet', 'CB+RFP', @(x) ischar(x) || isstring(x));
addParameter(parser, 'Tolerance', 1e-9, @(x) isnumeric(x) && isscalar(x) && x >= 0);
addParameter(parser, 'CacheFile', '', @(x) ischar(x) || isstring(x));
parse(parser, varargin{:});

sheetName = char(parser.Results.Sheet);
tolerance = parser.Results.Tolerance;
cacheFile = parser.Results.CacheFile;

data = import_CB_quantifications(workbookFile, 'Sheet', sheetName, ...
    'IncludeZeroCounts', true, 'CacheFile', cacheFile, 'ReimportExcel', true);
rawCells = readcell(char(data.workbookFile), 'Sheet', sheetName, 'Range', 'A:AQ');
sheetValues = cells_to_double_matrix(rawCells);

categories = strings(0, 1);
details = strings(0, 1);
mice = nan(0, 1);
workbookRows = nan(0, 1);
workbookColumns = strings(0, 1);
workbookValues = nan(0, 1);
expectedValues = nan(0, 1);
differences = nan(0, 1);
passed = false(0, 1);

rawTable = data.raw;
stackChecks = {
    'P', 'TotalRfp', 'total RFP'
    'Q', 'TotalOriens', 'oriens/above'
    'R', 'TotalDeep', 'deep'
    'S', 'TotalSuperficial', 'superficial'
    'T', 'TotalRadiatum', 'radiatum/below'
    'V', 'TotalRfp', 'total CB classes'
    'W', 'TotalCbPositive', 'CB+'
    'X', 'TotalCbNegative', 'CB-'
    'Y', 'TotalCbAmbiguous', 'CB+/-'};

for rowIdx = 1:height(rawTable)
    for checkIdx = 1:size(stackChecks, 1)
        col = stackChecks{checkIdx, 1};
        varName = stackChecks{checkIdx, 2};
        label = stackChecks{checkIdx, 3};
        add_check('stack total', label, rawTable.Mouse(rowIdx), rawTable.ExcelRow(rowIdx), ...
            col, value_at(sheetValues, rawTable.ExcelRow(rowIdx), col), rawTable.(varName)(rowIdx));
    end
end

summaryRows = find_mouse_summary_rows(sheetValues, rawTable);
summaryChecks = {
    'P', 'TotalRfp', 'total RFP'
    'Q', 'TotalOriens', 'oriens/above'
    'R', 'TotalDeep', 'deep'
    'S', 'TotalSuperficial', 'superficial'
    'T', 'TotalRadiatum', 'radiatum/below'
    'V', 'TotalRfp', 'total CB classes'
    'W', 'TotalCbPositive', 'CB+ total'
    'X', 'TotalCbNegative', 'CB- total'
    'Y', 'TotalCbAmbiguous', 'CB+/- total'
    'AA', 'TotalCbPositive', 'CB+ total'
    'AB', 'CbPositiveRadiatum', 'CB+ radiatum'
    'AC', 'CbPositiveSuperficial', 'CB+ superficial'
    'AD', 'CbPositiveDeep', 'CB+ deep'
    'AE', 'CbPositiveOriens', 'CB+ oriens'
    'AG', 'TotalCbAmbiguous', 'CB+/- total'
    'AH', 'CbAmbiguousRadiatum', 'CB+/- radiatum'
    'AI', 'CbAmbiguousSuperficial', 'CB+/- superficial'
    'AJ', 'CbAmbiguousDeep', 'CB+/- deep'
    'AK', 'CbAmbiguousOriens', 'CB+/- oriens'
    'AM', 'TotalCbNegative', 'CB- total'
    'AN', 'CbNegativeRadiatum', 'CB- radiatum'
    'AO', 'CbNegativeSuperficial', 'CB- superficial'
    'AP', 'CbNegativeDeep', 'CB- deep'
    'AQ', 'CbNegativeOriens', 'CB- oriens'};

for rowIdx = 1:numel(summaryRows)
    summaryRow = summaryRows(rowIdx);
    mouse = value_at(sheetValues, summaryRow, 'O');
    mouseMask = rawTable.Mouse == mouse;
    for checkIdx = 1:size(summaryChecks, 1)
        col = summaryChecks{checkIdx, 1};
        varName = summaryChecks{checkIdx, 2};
        label = summaryChecks{checkIdx, 3};
        expected = sum(rawTable.(varName)(mouseMask));
        add_check('mouse total', label, mouse, summaryRow, col, ...
            value_at(sheetValues, summaryRow, col), expected);
    end
end

pctTotalChecks = summaryChecks(~ismember(summaryChecks(:, 1), {'P', 'AA', 'AG', 'AM'}), :);
for rowIdx = 1:numel(summaryRows)
    summaryRow = summaryRows(rowIdx);
    pctRow = summaryRow + 7;
    mouse = value_at(sheetValues, summaryRow, 'O');
    mouseMask = rawTable.Mouse == mouse;
    totalRfp = sum(rawTable.TotalRfp(mouseMask));
    for checkIdx = 1:size(pctTotalChecks, 1)
        col = pctTotalChecks{checkIdx, 1};
        varName = pctTotalChecks{checkIdx, 2};
        label = pctTotalChecks{checkIdx, 3};
        expected = safe_percent(sum(rawTable.(varName)(mouseMask)), totalRfp);
        add_check('percent of total RFP', label, mouse, pctRow, col, ...
            value_at(sheetValues, pctRow, col), expected);
    end
end

classPercentChecks = {
    'AA', 'TotalCbPositive', 'TotalCbPositive', 'CB+ scale'
    'AB', 'CbPositiveRadiatum', 'TotalCbPositive', 'CB+ radiatum'
    'AC', 'CbPositiveSuperficial', 'TotalCbPositive', 'CB+ superficial'
    'AD', 'CbPositiveDeep', 'TotalCbPositive', 'CB+ deep'
    'AE', 'CbPositiveOriens', 'TotalCbPositive', 'CB+ oriens'
    'AG', 'TotalCbAmbiguous', 'TotalCbAmbiguous', 'CB+/- scale'
    'AH', 'CbAmbiguousRadiatum', 'TotalCbAmbiguous', 'CB+/- radiatum'
    'AI', 'CbAmbiguousSuperficial', 'TotalCbAmbiguous', 'CB+/- superficial'
    'AJ', 'CbAmbiguousDeep', 'TotalCbAmbiguous', 'CB+/- deep'
    'AK', 'CbAmbiguousOriens', 'TotalCbAmbiguous', 'CB+/- oriens'
    'AM', 'TotalCbNegative', 'TotalCbNegative', 'CB- scale'
    'AN', 'CbNegativeRadiatum', 'TotalCbNegative', 'CB- radiatum'
    'AO', 'CbNegativeSuperficial', 'TotalCbNegative', 'CB- superficial'
    'AP', 'CbNegativeDeep', 'TotalCbNegative', 'CB- deep'
    'AQ', 'CbNegativeOriens', 'TotalCbNegative', 'CB- oriens'};

scaleColumns = ["AA", "AG", "AM"];
for rowIdx = 1:numel(summaryRows)
    summaryRow = summaryRows(rowIdx);
    pctRow = summaryRow + 14;
    mouse = value_at(sheetValues, summaryRow, 'O');
    mouseMask = rawTable.Mouse == mouse;
    for checkIdx = 1:size(classPercentChecks, 1)
        col = classPercentChecks{checkIdx, 1};
        varName = classPercentChecks{checkIdx, 2};
        denomName = classPercentChecks{checkIdx, 3};
        label = classPercentChecks{checkIdx, 4};
        denom = sum(rawTable.(denomName)(mouseMask));
        if any(scaleColumns == string(col))
            expected = 100 / denom;
        else
            expected = safe_percent(sum(rawTable.(varName)(mouseMask)), denom);
        end
        add_check('percent within CB class', label, mouse, pctRow, col, ...
            value_at(sheetValues, pctRow, col), expected);
    end
end

report = table(categories, details, mice, workbookRows, workbookColumns, ...
    workbookValues, expectedValues, differences, passed, ...
    'VariableNames', {'Category', 'Detail', 'Mouse', 'WorkbookRow', ...
    'WorkbookColumn', 'WorkbookValue', 'ExpectedValue', 'Difference', 'Passed'});

if nargout == 0
    failed = report(~report.Passed, :);
    if isempty(failed)
        fprintf('All %d CB worksheet checks passed within tolerance %.3g.\n', height(report), tolerance);
    else
        fprintf('%d of %d CB worksheet checks failed within tolerance %.3g.\n', ...
            height(failed), height(report), tolerance);
        disp(failed(:, {'Category', 'Detail', 'Mouse', 'WorkbookRow', ...
            'WorkbookColumn', 'WorkbookValue', 'ExpectedValue', 'Difference'}));
    end
end

    function add_check(category, detail, mouse, workbookRow, workbookColumn, workbookValue, expectedValue)
        categories(end + 1, 1) = string(category);
        details(end + 1, 1) = string(detail);
        mice(end + 1, 1) = mouse;
        workbookRows(end + 1, 1) = workbookRow;
        workbookColumns(end + 1, 1) = string(workbookColumn);
        workbookValues(end + 1, 1) = workbookValue;
        expectedValues(end + 1, 1) = expectedValue;
        differences(end + 1, 1) = workbookValue - expectedValue;
        passed(end + 1, 1) = values_match(workbookValue, expectedValue, tolerance);
    end

end

function rows = find_mouse_summary_rows(sheetValues, rawTable)
rowNumbers = (1:size(sheetValues, 1))';
mouseColumn = sheetValues(:, excel_col_to_num('O'));
rows = rowNumbers(rowNumbers > max(rawTable.ExcelRow) & rowNumbers < 100 & ~isnan(mouseColumn));
end

function value = value_at(sheetValues, row, col)
value = sheetValues(row, excel_col_to_num(col));
end

function tf = values_match(a, b, tolerance)
tf = (isnan(a) && isnan(b)) || abs(a - b) <= tolerance;
end

function value = safe_percent(numerator, denominator)
if denominator == 0
    value = nan;
else
    value = numerator * 100 / denominator;
end
end

function matrix = cells_to_double_matrix(cells)
matrix = nan(size(cells));
for rowIdx = 1:size(cells, 1)
    for colIdx = 1:size(cells, 2)
        matrix(rowIdx, colIdx) = cell_to_double(cells{rowIdx, colIdx});
    end
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

function idx = excel_col_to_num(col)
col = char(upper(string(col)));
idx = 0;
for charIdx = 1:numel(col)
    idx = idx * 26 + double(col(charIdx)) - double('A') + 1;
end
end

function tf = is_name_value_start(value)
if ~(ischar(value) || (isstring(value) && isscalar(value)))
    tf = false;
    return
end
name = lower(char(string(value)));
validNames = {'sheet', 'tolerance', 'cachefile'};
tf = any(strcmp(name, validNames));
end
